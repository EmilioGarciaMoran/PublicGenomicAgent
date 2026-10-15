from __future__ import annotations

import json
from pathlib import Path

from .base import CaseManifest, GenomicRange, IndividualSpec


# ---------------------------------------------------------------------
# Carga y guardado
# ---------------------------------------------------------------------

def load_case_manifest(path: Path) -> CaseManifest:
    """Carga un CaseManifest desde JSON."""
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"CaseManifest no encontrado: {p}")
    data = json.loads(p.read_text())
    manifest = CaseManifest(**data)
    manifest.validate_minimum()
    return manifest


def save_case_manifest(manifest: CaseManifest, path: Path) -> Path:
    """Guarda un CaseManifest a JSON."""
    p = Path(path).expanduser().resolve()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(manifest.model_dump(), indent=2))
    return p


# ---------------------------------------------------------------------
# Conversores
# ---------------------------------------------------------------------

def case_manifest_to_pedigree(
    manifest: CaseManifest,
) -> list[IndividualSpec]:
    """Extrae el pedigrí del manifiesto.

    Si no hay pedigrí, lanza error. El pipeline debe verificar
    manifest.has_pedigree() antes de invocar tools que lo requieran.
    """
    if not manifest.has_pedigree():
        raise ValueError(
            f"CaseManifest '{manifest.case_id}' no tiene pedigrí completo "
            "(pedigree + proband). No se puede aplicar análisis mendeliano."
        )
    return list(manifest.pedigree)


def case_manifest_to_rois(manifest: CaseManifest) -> list[GenomicRange]:
    """Extrae los ROIs del manifiesto.

    Si no hay ROIs pero hay fenotipo, el sistema debería derivarlos
    consultando HPO -> genes. Esa lógica vive fuera de esta función.
    """
    if not manifest.has_rois():
        raise ValueError(
            f"CaseManifest '{manifest.case_id}' no tiene candidate_rois. "
            "El sistema debe derivarlos (por HPO o por búsqueda) antes "
            "de invocar tools que los requieran."
        )
    return list(manifest.candidate_rois)


def case_manifest_to_hpo_terms(manifest: CaseManifest) -> dict[str, list[str]]:
    """Devuelve el mapa sample -> [HPO IDs]."""
    if not manifest.has_phenotype():
        raise ValueError(
            f"CaseManifest '{manifest.case_id}' no tiene fenotipo. "
            "No se puede priorizar por similitud fenotípica."
        )
    return dict(manifest.hpo_terms)


# ---------------------------------------------------------------------
# Validación funcional — qué se puede y qué no con este manifiesto
# ---------------------------------------------------------------------

def capabilities(manifest: CaseManifest) -> dict[str, bool]:
    """Devuelve qué análisis se pueden hacer con este manifiesto.

    Útil para que el agente decida qué tools invocar.
    """
    return {
        "fetch_roi": manifest.has_rois(),
        "call_variants": manifest.has_rois(),      # necesita un ROI
        "mendelian_filter": manifest.has_pedigree(),
        "plink_validate": manifest.has_pedigree(),
        "phenotype_ranking": manifest.has_phenotype(),
        "consistency_check": (
            manifest.has_phenotype()
            and len(manifest.affected_samples) > 1
        ),
        "local_pangenome": manifest.has_rois(),
        "build_local_graph": manifest.has_rois(),
    }


# ---------------------------------------------------------------------
# Constructores desde otras fuentes
# ---------------------------------------------------------------------

def from_fam(
    fam_path: Path,
    case_id: str,
    proband: str,
    hpo_terms: dict[str, list[str]] | None = None,
    candidate_rois: list[GenomicRange] | None = None,
) -> CaseManifest:
    """Construye un CaseManifest a partir de un .fam de PLINK."""
    from .pedigree import load_pedigree_from_fam

    pedigree = load_pedigree_from_fam(fam_path)
    affected = [ind.sample for ind in pedigree if ind.affected]

    return CaseManifest(
        case_id=case_id,
        source="fam-file",
        pedigree=pedigree,
        proband=proband,
        hpo_terms=hpo_terms or {},
        affected_samples=affected,
        candidate_rois=candidate_rois or [],
        extractor="fam-file",
        confidence=1.0,
    )


def minimal_from_roi(
    case_id: str,
    chrom: str,
    start: int,
    end: int,
    label: str = "",
) -> CaseManifest:
    """Caso mínimo: solo un ROI, sin pedigrí ni fenotipo.

    Útil para estudios poblacionales, análisis exploratorios, o
    validación de variantes conocidas.
    """
    return CaseManifest(
        case_id=case_id,
        source="population",
        candidate_rois=[GenomicRange(
            chrom=chrom, start=start, end=end, label=label,
        )],
        extractor="manual",
        confidence=1.0,
    )


def minimal_from_phenotype(
    case_id: str,
    sample: str,
    hpo_terms: list[str],
    phenotype_text: str = "",
) -> CaseManifest:
    """Caso mínimo: solo fenotipo, sin pedigrí ni ROI.

    El sistema debe derivar ROIs por HPO -> genes.
    """
    return CaseManifest(
        case_id=case_id,
        source="bibliography",
        hpo_terms={sample: hpo_terms},
        phenotype_text={sample: phenotype_text} if phenotype_text else {},
        affected_samples=[sample],
        extractor="manual",
        confidence=1.0,
    )


# ---------------------------------------------------------------------
# Extracción desde lenguaje natural
# ---------------------------------------------------------------------

import re as _re


def _extract_json_block(raw: str) -> dict:
    """Extrae el primer bloque JSON válido del output del LLM.

    Robusto ante texto adicional, fences de markdown, o prefijos raros.
    """
    # Quitar fences de markdown ```json ... ```
    text = raw.strip()
    if text.startswith("```"):
        # Quitar la primera línea (```json o ```)
        lines = text.splitlines()
        text = "\n".join(lines[1:])
        if text.rstrip().endswith("```"):
            text = text.rstrip()[:-3]

    # Intentar parsear directo
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # Buscar el primer objeto {...} balanceado
    depth = 0
    start = None
    for i, ch in enumerate(text):
        if ch == "{":
            if depth == 0:
                start = i
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0 and start is not None:
                candidate = text[start : i + 1]
                try:
                    return json.loads(candidate)
                except json.JSONDecodeError:
                    start = None
                    depth = 0
    raise ValueError(
        f"No se encontró un objeto JSON válido en la respuesta del LLM. "
        f"Primeros 200 caracteres: {raw[:200]!r}"
    )


def _verify_samples_in_text(
    manifest: CaseManifest,
    text: str,
) -> tuple[list[str], list[str]]:
    """Verifica que los samples del manifiesto aparecen en el texto.

    Nivel 2 de validación: detecta si el LLM inventó samples.

    Returns:
        (verified, unverified): listas de samples.
    """
    text_lower = text.lower()
    verified: list[str] = []
    unverified: list[str] = []

    # Aceptamos "father" y "padre" como equivalentes, etc.
    SYNONYMS = {
        "father": ["father", "padre", "papá", "the father", "el padre"],
        "mother": ["mother", "madre", "mamá", "the mother", "la madre"],
        "proband": ["proband", "proband", "paciente", "patient", "hijo", "hija", "child"],
        "son": ["son", "hijo"],
        "daughter": ["daughter", "hija"],
        "sibling": ["sibling", "hermano", "hermana", "brother", "sister"],
    }

    for ind in manifest.pedigree:
        sample = ind.sample
        # Si el sample es literal en el texto, verificado
        if sample.lower() in text_lower:
            verified.append(sample)
            continue
        # Si el sample coincide con un sinónimo conocido, verificado
        synonyms = SYNONYMS.get(sample.lower(), [])
        if any(syn in text_lower for syn in synonyms):
            verified.append(sample)
            continue
        # Si no, no verificado
        unverified.append(sample)

    return verified, unverified


def _validate_and_clean(
    data: dict,
    text: str,
) -> tuple[CaseManifest, list[str]]:
    """Valida y limpia un dict como CaseManifest.

    Nivel 1: Pydantic valida el esquema.
    Nivel 2: verifica samples contra el texto.

    Returns:
        (manifest, warnings)
    """
    warnings: list[str] = []

    # Asegurar case_id
    if not data.get("case_id"):
        # Derivar uno simple
        import uuid
        data["case_id"] = f"nl_case_{uuid.uuid4().hex[:8]}"
        warnings.append(f"case_id derivado automáticamente: {data['case_id']}")

    # Asegurar extractor
    if not data.get("extractor"):
        data["extractor"] = "llm"

    # Limpiar campos que pueden venir con valores inesperados
    # (el LLM a veces pone null en lugar de lista vacía)
    for k in ("pedigree", "hpo_terms", "phenotype_text",
              "candidate_rois", "candidate_genes", "affected_samples"):
        if data.get(k) is None:
            data[k] = [] if k != "hpo_terms" and k != "phenotype_text" else {}

    # Pydantic valida
    manifest = CaseManifest.model_validate(data)

    # Nivel 2: verificar samples contra el texto
    verified, unverified = _verify_samples_in_text(manifest, text)
    if unverified:
        warnings.append(
            f"Los siguientes samples NO aparecen literalmente en el texto: "
            f"{unverified}. El LLM puede haberlos inventado."
        )

    return manifest, warnings


def nl_to_case_manifest(
    runtime: "ToolRuntime",
    inp: "NLToCaseManifestInput",
) -> "NLToCaseManifestOutput":
    """Extrae un CaseManifest de un relato clínico con un LLM local.

    No usa ToolRuntime para ejecutar binarios; lo recibe por consistencia
    con la firma de dispatch() pero no lo utiliza. La extracción la hace
    el LLM configurado (por defecto qwen2.5:7b en Ollama local).
    """
    from ..agent.config import load_config
    from ..agent.llm_factory import build_llm_client
    from ..agent.prompts import (
        EXTRACT_MANIFEST_SYSTEM_PROMPT,
        render_extract_manifest_user,
    )
    from .base import NLToCaseManifestInput, NLToCaseManifestOutput

    text = inp.text.strip()
    if not text:
        raise ValueError("El relato clínico está vacío")

    # Cliente LLM
    cfg = load_config()
    client = build_llm_client(cfg.llm)

    # Prompt
    user = render_extract_manifest_user(text, case_id=inp.case_id)

    # Llamada al LLM
    raw = client.complete(EXTRACT_MANIFEST_SYSTEM_PROMPT, user)

    # Parsear JSON
    data = _extract_json_block(raw)

    # Validar y limpiar
    manifest, warnings = _validate_and_clean(data, text)

    # Escribir a disco
    output_dir = Path(inp.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / f"{manifest.case_id}.manifest.json"
    save_case_manifest(manifest, manifest_path)

    # Verificar samples
    verified, unverified = _verify_samples_in_text(manifest, text)

    # Mensaje
    n_samples = len(manifest.pedigree)
    n_hpo = sum(len(v) for v in manifest.hpo_terms.values())
    n_rois = len(manifest.candidate_rois)
    n_genes = len(manifest.candidate_genes)
    message = (
        f"manifest '{manifest.case_id}': "
        f"{n_samples} sample(s), {n_hpo} HPO term(s), "
        f"{n_rois} ROI(s), {n_genes} gene(s)"
    )
    if warnings:
        message += f" ({len(warnings)} warning(s))"

    return NLToCaseManifestOutput(
        ok=True,
        tool="nl_to_case_manifest",
        message=message,
        outputs={"manifest": str(manifest_path)},
        manifest=manifest,
        manifest_json=manifest_path,
        samples_found=verified,
        samples_unverified=unverified,
        warnings=warnings,
    )
