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


def _key(k: Any, *, include_paths: bool = False) -> str:
    """Sanea una CLAVE del estado.

    Las claves del estado son etiquetas lógicas del tipo
    "called:father" o "roi_bam:proband". Si por algún motivo una clave
    contiene una ruta absoluta (p. ej. por un _register_side_effects
    antiguo), este helper la reduce a "prefijo:basename" antes de
    enviarla al LLM.

    Con `include_paths=True` se respeta la clave tal cual.
    """
    s = str(k)
    if include_paths:
        return s
    if "/" not in s:
        return s
    prefix, sep, rest = s.partition(":")
    if not sep:
        # Clave sin prefijo que contiene "/": reducir a basename
        return _basename(s)
    return f"{prefix}:{_basename(rest)}"


def render_state(state: AgentState, *, include_paths: bool = False) -> str:
    """Serializa el estado a un bloque de texto compacto para el LLM.

    Si `include_paths=False` (por defecto), las rutas se envían como
    basenames. Esto evita filtrar directorios absolutos con posible
    información identificativa del paciente.
    """
    case = state.case

    def _path(v: Any) -> str:
        return str(v) if include_paths else _basename(v)

    def _safe_key(k: Any) -> str:
        return _key(k, include_paths=include_paths)

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
            lines.append(f"  - {_safe_key(sample)}: {_path(p)}")
    if state.vcfs:
        lines.append("vcfs:")
        for label, p in state.vcfs.items():
            lines.append(f"  - {_safe_key(label)}: {_path(p)}")
    if state.references:
        lines.append("references:")
        for label, p in state.references.items():
            lines.append(f"  - {_safe_key(label)}: {_path(p)}")
    if state.artifacts:
        lines.append("artifacts:")
        for label, p in state.artifacts.items():
            lines.append(f"  - {_safe_key(label)}: {_path(p)}")

    # Historial de tools (resumido, sin argumentos completos)
    if state.tool_calls:
        lines.append("tool_calls_so_far:")
        for i, c in enumerate(state.tool_calls):
            status = "ok" if c.ok else f"error={c.error}"
            lines.append(f"  - [{i}] {c.tool_name}: {status}")

    return "\n".join(lines)


def render_tools_catalog(compact: bool = True) -> str:
    """Renderiza el catálogo de tools desde TOOL_REGISTRY.

    Por defecto, formato compacto: una línea por tool con su nombre,
    descripción y lista de campos required + opcionales. Esto reduce
    el prompt de ~7000 tokens a ~400, lo cual es crítico en MacAir
    con Qwen 7B (latencia y respeto de la instrucción de JSON puro).

    Con `compact=False`, vuelca el JSON Schema completo (útil solo
    para depurar).
    """
    if not compact:
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

    lines: list[str] = []
    for name in sorted(TOOL_REGISTRY):
        spec = TOOL_REGISTRY[name]
        schema = spec.input_schema()
        props = schema.get("properties", {})
        required = set(schema.get("required", []))
        fields = []
        for field_name, field_schema in props.items():
            ftype = field_schema.get("type", "any")
            if "anyOf" in field_schema:
                # Extraer el primer tipo no-null
                for alt in field_schema["anyOf"]:
                    if alt.get("type") != "null":
                        ftype = alt.get("type", "any")
                        break
            marker = "*" if field_name in required else ""
            fields.append(f"{field_name}{marker}:{ftype}")
        fields_str = ", ".join(fields)
        lines.append(f"- {name}: {spec.description} | fields: {fields_str}")
    return "\n".join(lines)


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


# ---------------------------------------------------------------------------
# Prompt para extraer CaseManifest de un relato clínico
# ---------------------------------------------------------------------------

EXTRACT_MANIFEST_SYSTEM_PROMPT = """\
Eres un asistente experto en genética clínica. Tu tarea es leer un
relato clínico y extraer la información estructurada que describe un
caso familiar, en formato JSON.

REGLAS DURAS (no negociables):

  1. Devuelve EXCLUSIVAMENTE un objeto JSON válido. Sin texto fuera
     del JSON. Sin markdown. Sin comentarios.

  2. NO INVENTES NADA. Si un dato no aparece explícitamente en el
     texto, déjalo vacío (lista vacía, dict vacío, null, o el default
     del campo). Es mejor un manifest incompleto y honesto que uno
     completo y falso.

  3. Los samples del pedigrí (father/mother/proband o nombres propios)
     deben aparecer LITERALMENTE en el texto. Si el texto dice "el
     padre" o "the father", usa el string "father" como sample. Si
     dice "Juan", usa "Juan".

  4. `consanguinity` es `true` solo si el texto menciona explícitamente
     parentesco entre los padres (primos, hermanos, tíos, "consanguíneo",
     "endogámico", "same tribe", "first cousins", etc.). En caso de duda,
     deja `false`.

  5. El campo `case_id` se deriva del texto si no se te da uno. Usa un
     slug corto y descriptivo (ej. "family_epilepsy_2026"). Sin espacios.

  6. Los HPO terms se expresan como IDs "HP:0000000". Si solo tienes el
     nombre del síntoma, NO inventes el ID: deja la lista vacía y añade
     el texto del síntoma a `phenotype_text`.

FORMATO DE SALIDA (JSON):

{
  "case_id": "string",
  "source": "clinical" | "bibliography" | "cohort" | "population",
  "pedigree": [
    {
      "sample": "string",
      "sex": "M" | "F" | "U",
      "affected": true | false,
      "father": "sample_or_null",
      "mother": "sample_or_null"
    }
  ],
  "proband": "sample_or_null",
  "consanguinity": true | false,
  "affected_samples": ["sample", ...],
  "hpo_terms": {"sample": ["HP:0000000", ...]},
  "phenotype_text": {"sample": "texto literal del relato para ese sample"},
  "candidate_rois": [],
  "candidate_genes": ["GEN1", ...],
  "confidence": 0.0-1.0,
  "extractor": "llm:qwen2.5:7b"
}

NOTAS:

  - `sex`: "M" (masculino), "F" (femenino), "U" (desconocido).
  - `father` / `mother`: el string del sample del padre/madre si están
    en el pedigrí, o `null`.
  - `confidence`: tu confianza global en la extracción. Si el texto es
    ambiguo o incompleto, pon un valor bajo (0.3-0.6).
  - El campo `source`:
      * "clinical" si es un caso clínico real.
      * "bibliography" si es de un paper.
      * "cohort" si menciona múltiples familias de un estudio.
      * "population" por defecto.
"""


def render_extract_manifest_user(
    text: str,
    case_id: str | None = None,
) -> str:
    """Renderiza el prompt de usuario para extracción de CaseManifest."""
    parts: list[str] = []
    if case_id:
        parts.append(f"case_id sugerido: {case_id}")
    parts.append("Relato clínico:")
    parts.append("---")
    parts.append(text.strip())
    parts.append("---")
    parts.append("Devuelve el JSON del CaseManifest. Nada más.")
    return "\n".join(parts)
