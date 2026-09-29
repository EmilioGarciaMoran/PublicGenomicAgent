"""Prompts y serialización de estado para el LLMPlanner.

Reglas de diseño:

  - El LLM NO recibe rutas absolutas por defecto. Solo basenames y
    metadatos. Esto es deliberado: en entornos clínicos (KSA, UE),
    las rutas absolutas pueden contener identificadores del paciente.
  - El catálogo de tools se renderiza desde TOOL_REGISTRY, así que
    siempre está sincronizado con el código. No hay prompt hardcodeado
    con una lista de tools que se pueda quedar obsoleta.
  - El SYSTEM_PROMPT es corto y con reglas duras. La creatividad del
    LLM no aporta aquí; lo que aporta es elegir la siguiente tool
    correcta dado el estado.

El contrato con el LLM es:

  ENTRADA:  system + user (user contiene estado + catálogo)
  SALIDA:   un único objeto JSON con la forma:
              {"tool_name": str, "args": {...}, "rationale": str}
            o bien
              {"stop": true, "rationale": str}
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..tools.registry import TOOL_REGISTRY
from .state import AgentState


SYSTEM_PROMPT = """\
Eres el planificador de PublicGenomicAgent, un agente de análisis
genómico clínico. Tu única tarea es decidir la SIGUIENTE acción a
ejecutar dado el estado actual de la sesión y el catálogo de tools
disponibles.

REGLAS DURAS (no negociables):

1. Solo puedes elegir tools que aparezcan en el catálogo.
2. Solo puedes proponer una tool si sus precondiciones se cumplen
   en el estado. Si faltan datos, responde con {"stop": true}.
3. Los argumentos deben ser coherentes con el schema de la tool.
   No inventes rutas ni valores: usa los que aparezcan en el estado.
4. No repitas una tool ya ejecutada con éxito para los mismos
   argumentos. Mira el historial de tool_calls.
5. Si ninguna tool aplica, responde {"stop": true, "rationale": "..."}.

FORMATO DE RESPUESTA: un único objeto JSON, sin markdown, sin
comentarios, sin texto adicional. Ejemplos válidos:

  {"tool_name": "qc_bam",
   "args": {"bam_path": "proband.bam", "level": "structural"},
   "rationale": "El BAM del probando aún no ha pasado QC."}

  {"stop": true, "rationale": "Todos los BAMs han pasado QC y no hay ROIs."}
"""


def _basename(p: Any) -> str:
    """Convierte una ruta a basename. Acepta Path o str."""
    if p is None:
        return ""
    try:
        return Path(str(p)).name
    except Exception:  # noqa: BLE001
        return str(p)


def render_state(state: AgentState, *, include_paths: bool = False) -> str:
    """Serializa el estado a un bloque de texto compacto para el LLM.

    Si `include_paths=False` (por defecto), las rutas se envían como
    basenames. Esto evita filtrar directorios absolutos con posible
    información identificativa del paciente.
    """
    case = state.case

    def _path(v: Any) -> str:
        return str(v) if include_paths else _basename(v)

    lines: list[str] = []
    lines.append(f"case_id: {case.case_id}")
    lines.append(f"source:  {case.source}")
    if case.proband:
        lines.append(f"proband: {case.proband}")
    lines.append(f"consanguinity: {case.consanguinity}")

    # Pedigrí
    if case.pedigree:
        lines.append("pedigree:")
        for ind in case.pedigree:
            rel = []
            if ind.father:
                rel.append(f"father={ind.father}")
            if ind.mother:
                rel.append(f"mother={ind.mother}")
            rel_s = (" " + " ".join(rel)) if rel else ""
            aff = "affected" if ind.affected else "unaffected"
            lines.append(f"  - {ind.sample} ({ind.sex}, {aff}){rel_s}")

    # ROIs
    if case.candidate_rois:
        lines.append("candidate_rois:")
        for roi in case.candidate_rois:
            label = f" {roi.label}" if roi.label else ""
            lines.append(f"  - {roi.chrom}:{roi.start}-{roi.end}{label}")

    # HPO
    if case.hpo_terms:
        lines.append("hpo_terms:")
        for sample, terms in case.hpo_terms.items():
            lines.append(f"  - {sample}: {', '.join(terms)}")

    # Recursos disponibles (basenames)
    if state.bams:
        lines.append("bams:")
        for sample, p in state.bams.items():
            lines.append(f"  - {sample}: {_path(p)}")
    if state.vcfs:
        lines.append("vcfs:")
        for label, p in state.vcfs.items():
            lines.append(f"  - {label}: {_path(p)}")
    if state.references:
        lines.append("references:")
        for label, p in state.references.items():
            lines.append(f"  - {label}: {_path(p)}")
    if state.artifacts:
        lines.append("artifacts:")
        for label, p in state.artifacts.items():
            lines.append(f"  - {label}: {_path(p)}")

    # Historial de tools (resumido, sin argumentos completos)
    if state.tool_calls:
        lines.append("tool_calls_so_far:")
        for i, c in enumerate(state.tool_calls):
            status = "ok" if c.ok else f"error={c.error}"
            lines.append(f"  - [{i}] {c.tool_name}: {status}")

    return "\n".join(lines)


def render_tools_catalog() -> str:
    """Renderiza el catálogo de tools desde TOOL_REGISTRY.

    Formato: JSON-lines, una tool por línea, con nombre, descripción,
    tags, si necesita runtime, y el JSON Schema del input.
    """
    out: list[str] = []
    for name in sorted(TOOL_REGISTRY):
        spec = TOOL_REGISTRY[name]
        entry = {
            "name": spec.name,
            "description": spec.description,
            "tags": list(spec.tags),
            "needs_runtime": spec.needs_runtime,
            "input_schema": spec.input_schema(),
        }
        out.append(json.dumps(entry, ensure_ascii=False))
    return "\n".join(out)


def render_user_prompt(
    state: AgentState, *, include_paths: bool = False
) -> str:
    """Construye el prompt de usuario: estado + catálogo."""
    return (
        "ESTADO ACTUAL:\n"
        f"{render_state(state, include_paths=include_paths)}\n\n"
        "CATÁLOGO DE TOOLS (JSON-lines):\n"
        f"{render_tools_catalog()}\n\n"
        "Devuelve la SIGUIENTE acción como un único objeto JSON."
    )
