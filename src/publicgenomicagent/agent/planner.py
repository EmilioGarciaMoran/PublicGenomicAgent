"""Planners: deciden la siguiente acción del agente.

Un `Planner` recibe un `AgentState` y devuelve un `PlannedAction`
(o `None` si no hay nada que hacer). El `AgentLoop` lo consulta en
cada iteración.

Dos implementaciones:

  - `RuleBasedPlanner`: reglas explícitas, sin LLM. Es el planner
    por defecto y también el fallback cuando el LLM falla.
  - `LLMPlanner`: usa un `LLMClient` inyectable. Valida la respuesta
    contra TOOL_REGISTRY y contra el input schema de la tool. Si algo
    falla, delega en un `RuleBasedPlanner` de fallback.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol

from ..tools.base import GenomicRange
from ..tools.registry import TOOL_REGISTRY, get_tool
from .llm import (
    LLMClient,
    LLMError,
    LLMResponse,
    LLMResponseFormatError,
    parse_llm_response,
)
from .prompts import SYSTEM_PROMPT, render_user_prompt
from .state import AgentState


# ---------------------------------------------------------------------------
# PlannedAction
# ---------------------------------------------------------------------------

@dataclass
class PlannedAction:
    """Acción propuesta por el planner: tool + args + justificación."""

    tool_name: str
    args: dict[str, Any]
    rationale: str = ""
    confidence: float = 1.0

    def describe(self) -> str:
        return f"{self.tool_name}({self.args}) — {self.rationale}"


class Planner(Protocol):
    """Contrato de un planner."""

    def next_action(self, state: AgentState) -> PlannedAction | None:
        """Devuelve la siguiente acción, o None si nada aplica."""
        ...


# ---------------------------------------------------------------------------
# RuleBasedPlanner
# ---------------------------------------------------------------------------

def _samples_in_case(state: AgentState) -> list[str]:
    if state.case.pedigree:
        return [ind.sample for ind in state.case.pedigree]
    return list(state.bams.keys())


def _roi_region(roi: GenomicRange) -> str:
    return f"{roi.chrom}:{roi.start}-{roi.end}"


def _roi_dir(state: AgentState) -> Path:
    d = state.artifacts.get("roi_dir")
    return Path(d) if d else Path.cwd() / "results" / "roi"


class RuleBasedPlanner:
    """Planner determinista basado en reglas explícitas.

    Prioridad:

      1. QC de cada BAM que aún no lo tenga.
      2. fetch_roi por (sample, ROI) que aún no exista.
      3. call_variants por ROI que aún no se haya llamado.
      4. mendelian_filter si hay pedigree + VCF de trío y no se ha hecho.

    Cada regla comprueba precondiciones en el estado antes de proponer
    la acción. Si ninguna aplica, devuelve None.
    """

    def __init__(self, roi_dir: Path | None = None):
        self._roi_dir_override = roi_dir

    def next_action(self, state: AgentState) -> PlannedAction | None:
        for rule in (
            self._rule_qc,
            self._rule_fetch_roi,
            self._rule_call_variants,
            self._rule_joint_call,
            self._rule_mendelian,
        ):
            action = rule(state)
            if action is not None:
                return action
        return None

    # --- reglas ---------------------------------------------------------

    def _rule_qc(self, state: AgentState) -> PlannedAction | None:
        for sample in _samples_in_case(state):
            bam = state.bams.get(sample)
            if bam is None:
                continue
            already = any(
                c.input.get("bam_path") == str(bam)
                and c.ok
                for c in state.calls_of("qc_bam")
            )
            if not already:
                return PlannedAction(
                    tool_name="qc_bam",
                    args={"bam_path": bam, "level": "structural"},
                    rationale=f"QC pendiente del BAM de {sample}.",
                )
        return None

    def _rule_fetch_roi(self, state: AgentState) -> PlannedAction | None:
        rois = state.case.candidate_rois
        samples = _samples_in_case(state)
        if not rois or not samples:
            return None

        done = {
            (c.input.get("bam_path"), c.input.get("region"))
            for c in state.calls_of("fetch_roi")
            if c.ok
        }
        out_dir = self._roi_dir_override or _roi_dir(state)

        for sample in samples:
            bam = state.bams.get(sample)
            if bam is None:
                continue
            for roi in rois:
                region = _roi_region(roi)
                if (str(bam), region) in done:
                    continue
                out_bam = out_dir / f"{sample}.{roi.chrom}_{roi.start}_{roi.end}.bam"
                return PlannedAction(
                    tool_name="fetch_roi",
                    args={
                        "bam_path": bam,
                        "region": region,
                        "output_bam": out_bam,
                    },
                    rationale=f"Falta sub-BAM de {sample} en {region}.",
                )
        return None

    def _rule_call_variants(self, state: AgentState) -> PlannedAction | None:
        rois = state.case.candidate_rois
        samples = _samples_in_case(state)
        if not rois or not samples:
            return None

        ref = state.references.get("hg38") or state.references.get("ref")
        if ref is None:
            return None

        out_dir = self._roi_dir_override or _roi_dir(state)
        called_pairs = {
            (c.input.get("region"), c.input.get("bam_path"))
            for c in state.calls_of("call_variants")
            if c.ok
        }

        for roi in rois:
            region = _roi_region(roi)
            # Buscamos el primer sample con sub-BAM ya disponible
            for sample in samples:
                sub_bam = out_dir / f"{sample}.{roi.chrom}_{roi.start}_{roi.end}.bam"
                if not sub_bam.exists():
                    continue
                if (region, str(sub_bam)) in called_pairs:
                    continue
                return PlannedAction(
                    tool_name="call_variants",
                    args={
                        "bam_path": sub_bam,
                        "reference_fasta": ref,
                        "output_vcf": sub_bam.with_suffix(".vcf.gz"),
                        "region": region,
                    },
                    rationale=f"Variant calling pendiente en {region} ({sample}).",
                )
        return None

    def _rule_joint_call(self, state: AgentState) -> PlannedAction | None:
        """Ejecuta `joint_call` una vez cuando los sub-BAMs del trío existen.

        Precondiciones:
          - Hay pedigree (caso familiar) y al menos 2 muestras.
          - Existe un sub-BAM por muestra para la primera ROI.
          - NO existe ya un VCF conjunto en state.vcfs["trio"].
          - NO se ha ejecutado joint_call antes para esa ROI.

        Si se cumplen, propone joint_call con todos los sub-BAMs del
        ROI en una sola llamada. Eso produce un VCF multi-sample con
        las columnas de muestra, que es lo que `mendelian_filter`
        necesita para verificar segregación.
        """
        if not state.case.has_pedigree():
            return None

        # ¿Ya hay VCF conjunto?
        if "trio" in state.vcfs:
            return None

        rois = state.case.candidate_rois
        samples = _samples_in_case(state)
        if not rois or len(samples) < 2:
            return None

        ref = state.references.get("hg38") or state.references.get("ref")
        if ref is None:
            return None

        out_dir = self._roi_dir_override or _roi_dir(state)

        # Buscar la primera ROI con sub-BAMs para todos los samples
        for roi in rois:
            region = _roi_region(roi)
            sub_bams: list[tuple[str, Path]] = []
            for sample in samples:
                sub_bam = out_dir / f"{sample}.{roi.chrom}_{roi.start}_{roi.end}.bam"
                if not sub_bam.exists():
                    break
                sub_bams.append((sample, sub_bam))
            else:
                # Todos los samples tienen sub-BAM para esta ROI.
                # ¿Se ha ejecutado joint_call para esta ROI?
                joint_calls = [
                    c for c in state.calls_of("joint_call") if c.ok
                ]
                region_done = any(
                    c.input.get("region") == region for c in joint_calls
                )
                if region_done:
                    continue

                out_vcf = out_dir / f"joint_{roi.chrom}_{roi.start}_{roi.end}.vcf.gz"
                return PlannedAction(
                    tool_name="joint_call",
                    args={
                        "bams": sub_bams,
                        "reference_fasta": ref,
                        "output_vcf": out_vcf,
                        "region": region,
                    },
                    rationale=(
                        f"Hay {len(sub_bams)} sub-BAMs del trío en {region}; "
                        f"generamos VCF conjunto para filtrar por mendelismo."
                    ),
                )
        return None

    def _rule_mendelian(self, state: AgentState) -> PlannedAction | None:
        if not state.case.has_pedigree():
            return None
        if state.has_tool("mendelian_filter"):
            return None
        trio_vcf = state.vcfs.get("trio") or state.vcfs.get("joint")
        if trio_vcf is None:
            return None
        return PlannedAction(
            tool_name="mendelian_filter",
            args={
                "trio_vcf": trio_vcf,
                "output_dir": Path.cwd() / "results" / "mendelian",
                "pedigree": [ind.model_dump() for ind in state.case.pedigree],
                "proband": state.case.proband,
            },
            rationale="Hay pedigree y VCF de trío, aplicamos filtros mendelianos.",
        )


# ---------------------------------------------------------------------------
# LLMPlanner
# ---------------------------------------------------------------------------

@dataclass
class LLMPlanner:
    """Planner basado en LLM con validación estricta y fallback.

    Flujo:

      1. Construye el prompt con `render_user_prompt(state)`.
      2. Llama a `client.complete(SYSTEM_PROMPT, user)`.
      3. Parsea el JSON con `parse_llm_response`.
      4. Valida el tool_name contra TOOL_REGISTRY.
      5. Valida los args contra el modelo Pydantic de la tool
         (con `dispatch(..., validate_only=True)`).
      6. Si algo falla en 2–5, delega en `fallback.next_action(state)`.
         Nunca reintenta ciegamente.

    El planner anota en `state.notes` cada decisión (éxito o fallback),
    para que la traza quede auditable.
    """

    client: LLMClient
    fallback: RuleBasedPlanner = field(default_factory=RuleBasedPlanner)
    include_paths: bool = False

    def next_action(self, state: AgentState) -> PlannedAction | None:
        user_prompt = render_user_prompt(state, include_paths=self.include_paths)

        try:
            raw = self.client.complete(SYSTEM_PROMPT, user_prompt)
        except LLMError as e:
            state.note(f"llm_planner: error del cliente ({e}); fallback")
            return self.fallback.next_action(state)
        except Exception as e:  # noqa: BLE001
            state.note(
                f"llm_planner: excepción inesperada ({type(e).__name__}: {e}); fallback"
            )
            return self.fallback.next_action(state)

        try:
            resp = parse_llm_response(raw)
        except LLMResponseFormatError as e:
            state.note(f"llm_planner: formato inválido ({e}); fallback")
            return self.fallback.next_action(state)

        if resp.stop:
            state.note(f"llm_planner: stop ({resp.rationale})")
            return None

        # Validación 1: tool conocida
        if resp.tool_name not in TOOL_REGISTRY:
            state.note(
                f"llm_planner: tool desconocida '{resp.tool_name}'; fallback"
            )
            return self.fallback.next_action(state)

        # Validación 2: args válidos según Pydantic
        spec = get_tool(resp.tool_name)
        try:
            spec.input_model(**resp.args)
        except Exception as e:  # noqa: BLE001
            state.note(
                f"llm_planner: args inválidos para {resp.tool_name} "
                f"({type(e).__name__}: {e}); fallback"
            )
            return self.fallback.next_action(state)

        state.note(
            f"llm_planner: {resp.tool_name} — {resp.rationale}"
        )
        return PlannedAction(
            tool_name=resp.tool_name,
            args=resp.args,
            rationale=resp.rationale,
            confidence=1.0,
        )

    # --- helpers ---------------------------------------------------------

    @staticmethod
    def _already_executed(
        state: AgentState,
        tool_name: str,
        args: dict,
    ) -> bool:
        """True si (tool_name, args) coincide con una tool_call ok previa.

        Corta el bucle infinito cuando el LLM insiste en la misma
        tool con los mismos argumentos.
        """
        def norm(d: dict) -> tuple:
            return tuple(sorted((k, str(v)) for k, v in d.items()))

        target_args = norm(args)

        for c in state.calls_of(tool_name):
            if not c.ok:
                continue
            if norm(c.input) == target_args:
                return True
        return False

# ---------------------------------------------------------------------------
# RuleFirstLLMPlanner (híbrido)
# ---------------------------------------------------------------------------

@dataclass
class RuleFirstLLMPlanner:
    """Planner híbrido: reglas primero, LLM como desambiguador.

    Diferencias con `LLMPlanner`:

      - `LLMPlanner` consulta SIEMPRE al LLM primero; solo cae al
        fallback si el LLM falla o alucina.
      - `RuleFirstLLMPlanner` consulta SIEMPRE al `RuleBasedPlanner`
        primero; solo llama al LLM si el determinismo no sabe qué
        hacer (devuelve None).

    Motivo del cambio: los modelos pequeños (3B) no siguen de forma
    fiable instrucciones complejas del tipo "no repitas tools ya
    ejecutadas" o "elige la siguiente tool según el historial".
    En pipelines con un orden canónico claro, el determinismo es
    superior; el LLM solo aporta valor cuando el estado es ambiguo.

    Casos en los que el LLM se consulta:
      - El RuleBasedPlanner no tiene regla aplicable (estado ambiguo).
      - Hay que desambiguar entre varias acciones válidas.
      - El caso es atípico y las reglas no lo cubren.

    Cuando el LLM se consulta, el flujo de validación es idéntico al
    del `LLMPlanner` (parseo, validación de tool, validación de args,
    rechazo de duplicados, resolución de paths).
    """

    client: LLMClient
    rule: RuleBasedPlanner = field(default_factory=RuleBasedPlanner)
    include_paths: bool = False

    def next_action(self, state: AgentState) -> PlannedAction | None:
        # 1. Intentar primero con reglas deterministas.
        rule_action = self.rule.next_action(state)
        if rule_action is not None:
            state.note(
                f"hybrid_planner: rule-based eligió "
                f"{rule_action.tool_name}; LLM no consultado"
            )
            return rule_action

        # 2. Reglas agotadas. Pedimos al LLM que desambigüe.
        state.note("hybrid_planner: reglas agotadas; consultando LLM")

        user_prompt = render_user_prompt(state, include_paths=self.include_paths)

        try:
            raw = self.client.complete(SYSTEM_PROMPT, user_prompt)
        except LLMError as e:
            state.note(f"hybrid_planner: error del cliente ({e}); stop")
            return None
        except Exception as e:  # noqa: BLE001
            state.note(
                f"hybrid_planner: excepción ({type(e).__name__}: {e}); stop"
            )
            return None

        try:
            resp = parse_llm_response(raw)
        except LLMResponseFormatError as e:
            state.note(f"hybrid_planner: formato inválido ({e}); stop")
            return None

        if resp.stop:
            state.note(f"hybrid_planner: stop ({resp.rationale})")
            return None

        if resp.tool_name not in TOOL_REGISTRY:
            state.note(
                f"hybrid_planner: tool desconocida '{resp.tool_name}'; stop"
            )
            return None

        spec = get_tool(resp.tool_name)
        try:
            spec.input_model(**resp.args)
        except Exception as e:  # noqa: BLE001
            state.note(
                f"hybrid_planner: args inválidos para {resp.tool_name} "
                f"({type(e).__name__}: {e}); stop"
            )
            return None

        # Resolver basenames a paths absolutos, luego comprobar
        # duplicados con los args ya resueltos.
        resolved_args = self._resolve_paths(state, resp.args)

        if self._already_executed(state, resp.tool_name, resolved_args):
            state.note(
                f"hybrid_planner: {resp.tool_name} ya ejecutada; stop"
            )
            return None

        state.note(f"hybrid_planner: LLM eligió {resp.tool_name} — {resp.rationale}")
        return PlannedAction(
            tool_name=resp.tool_name,
            args=resolved_args,
            rationale=resp.rationale,
            confidence=0.9,
        )

    # --- helpers (duplicados del LLMPlanner; si crece la duplicación,
    #     extraer a mixin) ------------------------------------------------

    @staticmethod
    def _already_executed(
        state: AgentState,
        tool_name: str,
        args: dict,
    ) -> bool:
        def norm(d: dict) -> tuple:
            return tuple(sorted((k, str(v)) for k, v in d.items()))

        target_args = norm(args)

        for c in state.calls_of(tool_name):
            if norm(c.input) == target_args:
                return True
        return False

    @staticmethod
    def _resolve_paths(state: AgentState, args: dict) -> dict:
        """Traduce basenames a paths absolutos usando el estado."""
        index: dict[str, Path] = {}

        def add(p) -> None:
            try:
                bp = Path(p)
                index[bp.name] = bp
            except Exception:  # noqa: BLE001
                pass

        for p in state.bams.values():
            add(p)
        for p in state.vcfs.values():
            add(p)
        for p in state.references.values():
            add(p)
        for p in state.artifacts.values():
            add(p)

        resolved: dict = {}
        for k, v in args.items():
            if isinstance(v, str) and v in index:
                resolved[k] = index[v]
            else:
                resolved[k] = v
        return resolved

