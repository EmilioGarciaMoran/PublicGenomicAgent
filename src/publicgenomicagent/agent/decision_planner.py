"""DecisionPlanner: planificación determinista con elección acotada del LLM.

Principio: el código decide QUÉ acciones son válidas (y construye sus
argumentos); el LLM, como mucho, elige ENTRE acciones ya válidas, y solo
cuando de verdad hay una bifurcación (pasos opcionales).

  - Sin candidatos            -> stop. El LLM NO se consulta.
  - Candidato requerido       -> se ejecuta. El LLM NO se consulta.
  - Solo candidatos opcionales-> bifurcación real: el LLM elige entre
                                 {candidatos..., "SKIP"} con salida
                                 restringida por JSON schema (enum).
  - Cualquier fallo del LLM   -> acción por defecto determinista
                                 ("fail closed"), anotada en el log.

El LLM nunca genera argumentos ni nombres de tool libres: solo devuelve
un valor del enum. Eso elimina la alucinación de tools/args/paths.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from typing import Any, Callable, Protocol

from .planner import PlannedAction, RuleBasedPlanner
from .state import AgentState

SKIP = "SKIP"

# Etapa canónica del pipeline y si el paso es obligatorio.
# (tool_name, stage, required)
_STAGES: tuple[tuple[str, int, bool], ...] = (
    ("qc_bam", 0, True),
    ("fetch_roi", 1, True),
    ("call_variants", 2, True),
    ("joint_call", 3, True),
    ("annotate_variants", 4, False),
    ("mendelian_filter", 5, True),
    ("extract_hpo", 6, False),
    ("phenotype_ranking", 7, False),
    ("prioritize_variants", 8, False),
)
_STAGE = {t: s for t, s, _ in _STAGES}
_REQUIRED = {t for t, _, r in _STAGES if r}

# regla del RuleBasedPlanner que genera la acción de cada tool
_RULE_ATTR = {
    "qc_bam": "_rule_qc",
    "fetch_roi": "_rule_fetch_roi",
    "call_variants": "_rule_call_variants",
    "joint_call": "_rule_joint_call",
    "annotate_variants": "_rule_annotate",
    "mendelian_filter": "_rule_mendelian",
    "extract_hpo": "_rule_extract_hpo",
    "phenotype_ranking": "_rule_phenotype_ranking",
    "prioritize_variants": "_rule_prioritize",
}

FORK_SYSTEM_PROMPT = (
    "You are a routing function inside a clinical genomics pipeline. "
    "Pick exactly one option from the allowed list. "
    "Prefer the first option unless the case summary clearly says it is "
    "not useful. Output JSON only."
)


class Chooser(Protocol):
    """Devuelve JSON con {"choice": <valor del enum>} para un schema dado."""

    model_id: str

    def choose(self, system: str, user: str, schema: dict) -> str: ...


@dataclass
class OllamaChooser:
    """Chooser con salida estructurada de Ollama (`format` = JSON schema).

    Usa `/api/generate`, igual que `OllamaClient`, y recupera el digest
    del modelo (best effort) para dejarlo en la traza de auditoría.
    """

    model: str = "qwen2.5:3b"
    endpoint: str = "http://localhost:11434/api/generate"
    timeout: float = 120.0
    num_ctx: int = 2048
    seed: int = 0
    _digest: str | None = field(default=None, init=False, repr=False)

    @property
    def model_id(self) -> str:
        return f"{self.model}@{self.model_digest() or 'unknown'}"

    def model_digest(self) -> str | None:
        if self._digest is None:
            import urllib.request

            base = self.endpoint.split("/api/")[0]
            try:
                with urllib.request.urlopen(f"{base}/api/tags", timeout=5) as r:
                    for m in json.loads(r.read()).get("models", []):
                        if m.get("name") == self.model or m.get("model") == self.model:
                            self._digest = str(m.get("digest", ""))[:12] or None
            except Exception:  # noqa: BLE001
                self._digest = None
        return self._digest

    def choose(self, system: str, user: str, schema: dict) -> str:
        import urllib.request

        body = {
            "model": self.model,
            "system": system,
            "prompt": user,
            "stream": False,
            "format": schema,
            "options": {"temperature": 0, "seed": self.seed,
                        "num_ctx": self.num_ctx, "num_predict": 64},
        }
        req = urllib.request.Request(
            self.endpoint,
            data=json.dumps(body).encode(),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=self.timeout) as r:
            return json.loads(r.read())["response"]


def build_chooser(cfg) -> "OllamaChooser | None":
    """Chooser a partir de `LLMConfig`. None si provider != ollama."""
    if (cfg.provider or "none").lower() != "ollama":
        return None
    return OllamaChooser(model=cfg.model, endpoint=cfg.endpoint,
                         timeout=cfg.timeout_seconds, num_ctx=cfg.num_ctx)


def _norm(d: dict) -> tuple:
    return tuple(sorted((k, str(v)) for k, v in d.items()))


@dataclass
class DecisionPlanner:
    """Planner determinista con bifurcaciones opcionales resueltas por LLM."""

    chooser: Chooser | None = None
    rules: RuleBasedPlanner = field(default_factory=RuleBasedPlanner)
    max_attempts: int = 2          # reintentos por (tool, args) fallidos
    decisions: list[dict[str, Any]] = field(default_factory=list)
    _skipped: set[str] = field(default_factory=set)

    # --- API Planner -----------------------------------------------------

    def next_action(self, state: AgentState) -> PlannedAction | None:
        cands = self._candidates(state)
        if not cands:
            return self._log(state, [], None, "rule", "sin candidatos: stop")

        first = cands[0]
        if first.tool_name in _REQUIRED:
            return self._log(state, cands, first, "rule", "paso requerido")

        # Bifurcación: opcionales por delante del siguiente requerido.
        next_req = next(
            (_STAGE[c.tool_name] for c in cands if c.tool_name in _REQUIRED),
            99,
        )
        fork = [c for c in cands
                if c.tool_name not in _REQUIRED and _STAGE[c.tool_name] < next_req]
        default = fork[0]

        if self.chooser is None or len(fork) == 0:
            return self._log(state, fork, default, "rule", "sin chooser: defecto")

        choice, src, note = self._ask(state, fork)
        if choice == SKIP:
            for c in fork:
                self._skipped.add(c.tool_name)
            return self.next_action(state)  # re-evalúa sin los saltados
        picked = next(c for c in fork if c.tool_name == choice)
        return self._log(state, fork, picked, src, note)

    # --- candidatos ------------------------------------------------------

    def _candidates(self, state: AgentState) -> list[PlannedAction]:
        out: list[PlannedAction] = []
        for tool, _stage, _req in _STAGES:
            if tool in self._skipped:
                continue
            action = getattr(self.rules, _RULE_ATTR[tool])(state)
            if action is None:
                continue
            if self._failures(state, action) >= self.max_attempts:
                state.note(f"decision_planner: {tool} agotó reintentos; descartado")
                continue
            out.append(action)
        return out

    @staticmethod
    def _failures(state: AgentState, action: PlannedAction) -> int:
        target = _norm(action.args)
        return sum(
            1 for c in state.calls_of(action.tool_name)
            if not c.ok and _norm(c.input) == target
        )

    # --- LLM restringido -------------------------------------------------

    def _ask(self, state: AgentState, fork: list[PlannedAction]):
        names = [c.tool_name for c in fork]
        enum = names + [SKIP]
        schema = {
            "type": "object",
            "properties": {"choice": {"type": "string", "enum": enum}},
            "required": ["choice"],
            "additionalProperties": False,
        }
        user = (
            "Case summary:\n" + self._summary(state)
            + "\nAllowed options:\n"
            + "\n".join(f"- {c.tool_name}: {c.rationale}" for c in fork)
            + f"\n- {SKIP}: skip these optional steps\n"
        )
        meta = {
            "model": getattr(self.chooser, "model_id", "?"),
            "prompt_sha256": hashlib.sha256(
                (FORK_SYSTEM_PROMPT + user).encode()).hexdigest()[:16],
        }
        try:
            raw = self.chooser.choose(FORK_SYSTEM_PROMPT, user, schema)
            choice = json.loads(raw)["choice"]
            if choice not in enum:
                raise ValueError(f"valor fuera del enum: {choice!r}")
        except Exception as e:  # noqa: BLE001  fail closed
            state.note(f"decision_planner: LLM falló ({type(e).__name__}); defecto")
            self.decisions.append({"event": "llm_fallback", "error": str(e), **meta})
            return names[0], "fallback", "defecto tras fallo del LLM"
        self.decisions.append({"event": "llm_choice", "choice": choice,
                               "options": enum, **meta})
        return choice, "llm", "elección acotada"

    @staticmethod
    def _summary(state: AgentState) -> str:
        done = sorted(t for t in _STAGE
                      if any(c.ok for c in state.calls_of(t)))
        return (f"done_tools={done}; pedigree={state.case.has_pedigree()}; "
                f"has_clinical_text={bool(state.case.phenotype_text)}")

    # --- log -------------------------------------------------------------

    def _log(self, state, cands, chosen, source, why):
        state.note(
            f"decision_planner[{source}]: "
            f"{chosen.tool_name if chosen else 'STOP'} ({why}); "
            f"candidatos={[c.tool_name for c in cands]}"
        )
        self.decisions.append({
            "event": "decision", "source": source,
            "chosen": chosen.tool_name if chosen else None,
            "candidates": [c.tool_name for c in cands],
        })
        return chosen
