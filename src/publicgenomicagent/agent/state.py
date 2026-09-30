"""Estado tipado y auditable de una sesión del agente.

Un `AgentState` representa todo lo que el agente sabe en un momento
dado: el manifiesto del caso, las rutas de BAMs y VCFs disponibles,
los artefactos generados por tools previas, y un log completo de cada
tool invocada.

`SessionContext` envuelve `AgentState` + `ToolRuntime` y expone
`call_tool(name, **kwargs)` como única puerta de entrada para ejecutar
tools. Esto garantiza que toda llamada quede registrada con su input
validado, su output serializado, su duración y, si falla, el error.
"""
from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from ..env.runtime import ToolRuntime
from ..tools.base import CaseManifest, ToolOutput
from ..tools.registry import dispatch, get_tool


class ToolCallRecord(BaseModel):
    """Traza auditable de una llamada a una tool."""

    tool_name: str
    input: dict[str, Any]
    output: dict[str, Any] | None = None
    ok: bool = True
    duration_ms: int = 0
    timestamp: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    error: str | None = None


class AgentState(BaseModel):
    """Estado completo de una sesión del agente.

    Mutable por diseño: el loop agéntico irá añadiendo BAMs, VCFs,
    artefactos y registros de tools conforme avanza el análisis.
    """

    case: CaseManifest

    # Rutas conocidas de la sesión, indexadas por nombre lógico.
    bams: dict[str, Path] = Field(default_factory=dict)     # sample -> BAM
    vcfs: dict[str, Path] = Field(default_factory=dict)     # sample/label -> VCF
    fastqs: dict[str, Path] = Field(default_factory=dict)   # sample -> FASTQ
    references: dict[str, Path] = Field(default_factory=dict)  # label -> FASTA

    # Artefactos intermedios generados por tools (grafos vg, packs, índices).
    artifacts: dict[str, Path] = Field(default_factory=dict)

    # Trazas y resultados.
    tool_calls: list[ToolCallRecord] = Field(default_factory=list)
    outputs: dict[str, dict[str, Any]] = Field(default_factory=dict)
    notes: list[str] = Field(default_factory=list)

    # --- API de consulta -------------------------------------------------

    def latest(self, tool_name: str) -> dict[str, Any] | None:
        """Último output registrado para una tool, o None."""
        return self.outputs.get(tool_name)

    def calls_of(self, tool_name: str) -> list[ToolCallRecord]:
        return [c for c in self.tool_calls if c.tool_name == tool_name]

    def has_tool(self, tool_name: str) -> bool:
        return tool_name in self.outputs

    def note(self, message: str) -> None:
        self.notes.append(message)

    # --- API de registro -------------------------------------------------

    def record(
        self,
        *,
        tool_name: str,
        inp: dict[str, Any],
        out: dict[str, Any] | None,
        duration_ms: int,
        error: str | None = None,
    ) -> ToolCallRecord:
        ok = error is None
        rec = ToolCallRecord(
            tool_name=tool_name,
            input=inp,
            output=out,
            ok=ok,
            duration_ms=duration_ms,
            error=error,
        )
        self.tool_calls.append(rec)
        if ok and out is not None:
            self.outputs[tool_name] = out
        return rec

    # --- Persistencia ----------------------------------------------------

    def to_json(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(
                self.model_dump(mode="json"),
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
            )
        )
        return path


class SessionContext:
    """Contexto vivo de una sesión: estado + runtime.

    Uso típico desde el loop agéntico:

        ctx = SessionContext(state=state, runtime=runtime)
        out = ctx.call_tool("fetch_roi",
                            bam_path=..., region=..., output_bam=...)
        # el estado queda con la traza actualizada
    """

    def __init__(self, state: AgentState, runtime: ToolRuntime | None = None):
        self.state = state
        self.runtime = runtime

    def call_tool(self, tool_name: str, **kwargs: Any) -> ToolOutput:
        """Ejecuta una tool por nombre, registra la traza y devuelve el output.

        - Valida el input contra el modelo Pydantic de la tool.
        - Si la tool necesita runtime y no lo hay, lanza ValueError.
        - Cualquier excepción de la tool se registra como error y se re-lanza.
        """
        spec = get_tool(tool_name)

        # Validación previa para poder registrar el input ya normalizado
        # (Path -> str, Enums -> valor, etc.) aunque la ejecución falle.
        try:
            validated = dispatch(tool_name, validate_only=True, **kwargs)
            inp_json = validated.model_dump(mode="json")
        except Exception as e:  # noqa: BLE001
            # Input inválido: lo registramos con error y re-lanzamos.
            self.state.record(
                tool_name=tool_name,
                inp={k: _jsonable(v) for k, v in kwargs.items()},
                out=None,
                duration_ms=0,
                error=f"input validation failed: {type(e).__name__}: {e}",
            )
            raise

        t0 = time.perf_counter()
        try:
            out = dispatch(tool_name, runtime=self.runtime, **kwargs)
        except Exception as e:  # noqa: BLE001
            elapsed = int((time.perf_counter() - t0) * 1000)
            self.state.record(
                tool_name=tool_name,
                inp=inp_json,
                out=None,
                duration_ms=elapsed,
                error=f"{type(e).__name__}: {e}",
            )
            raise

        elapsed = int((time.perf_counter() - t0) * 1000)
        out_json = out.model_dump(mode="json") if hasattr(out, "model_dump") else None
        self.state.record(
            tool_name=tool_name,
            inp=inp_json,
            out=out_json,
            duration_ms=elapsed,
        )

        # Efectos colaterales útiles: si la tool produce ficheros
        # convencionales, los indexamos en el estado.
        _register_side_effects(self.state, tool_name, out)

        return out


def _jsonable(value: Any) -> Any:
    """Convierte Path, Enum, etc. a algo serializable a JSON."""
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "value"):  # Enum
        return value.value
    return value


def _logical_label(path: str | Path) -> str:
    """Extrae una etiqueta lógica legible de un path.

    Ejemplos:
      /r/roi/father.chr1_1000_1200.vcf.gz  -> father
      /r/roi/proband.chr1_1000_1200.bam    -> proband
      /r/graph.vg                          -> graph
      /r/pack.pack                         -> pack

    Estrategia:
      1. Tomar el basename.
      2. Quitar la extensión compuesta (.vcf.gz, .bam.bai, .bam, .vg, ...).
      3. Cortar por "." y devolver el primer trozo.

    Nota: el corte por "." asume que el primer componente del nombre
    es el identificador lógico (father, proband, graph, pack). Si en
    el futuro un nombre de fichero legítimamente empieza por puntos
    (p. ej. "sample.s1.ROI.bam" donde sample.s1 es la etiqueta),
    habría que ajustar esta función. Hoy por hoy el prefijo es
    siempre el sample ID o el tipo de artefacto.
    """
    name = Path(str(path)).name
    # Quitar extensiones compuestas y simples.
    for suffix in (
        ".vcf.gz.tbi",
        ".vcf.gz",
        ".bam.bai",
        ".cram.crai",
        ".bam",
        ".cram",
        ".vcf",
        ".bcf",
        ".vg",
        ".xg",
        ".pack",
        ".gam",
        ".fa",
        ".fai",
    ):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    else:
        # Sin extensión conocida: quitar solo el último componente
        if "." in name:
            name = name.rsplit(".", 1)[0]
    # Primer trozo del nombre separado por "."
    return name.split(".")[0]


def _register_side_effects(state: AgentState, tool_name: str, out: ToolOutput) -> None:
    """Indexa rutas producidas por tools en el estado de la sesión.

    Las CLAVES son etiquetas lógicas (p. ej. "called:father"), no rutas
    absolutas. Esto es deliberado: las claves se envían al LLM en
    `render_state`, y no deben filtrar directorios clínicos.

    Los VALORES sí son Path absolutos, porque el propio agente los
    necesita para ejecutar tools. La sanitización de valores ocurre
    en `render_state` (a basename) salvo que include_paths=True.
    """
    d = out.model_dump(mode="json") if hasattr(out, "model_dump") else {}

    if tool_name == "fetch_roi":
        if "output_bam" in d:
            label = _logical_label(d["output_bam"])
            state.artifacts[f"roi_bam:{label}"] = Path(d["output_bam"])
    elif tool_name == "call_variants":
        if "output_vcf" in d:
            label = _logical_label(d["output_vcf"])
            state.vcfs[f"called:{label}"] = Path(d["output_vcf"])
    elif tool_name == "build_local_graph":
        if "graph_vg" in d:
            state.artifacts["graph_vg"] = Path(d["graph_vg"])
        if "graph_xg" in d:
            state.artifacts["graph_xg"] = Path(d["graph_xg"])
    elif tool_name == "align_to_graph":
        if d.get("pack"):
            state.artifacts["pack"] = Path(d["pack"])
    elif tool_name == "call_from_graph":
        if "vcf" in d:
            label = _logical_label(d["vcf"])
            state.vcfs[f"graph_called:{label}"] = Path(d["vcf"])
    elif tool_name == "local_pangenome":
        if "enriched_fasta" in d:
            label = _logical_label(d["enriched_fasta"])
            state.references[f"enriched:{label}"] = Path(d["enriched_fasta"])
