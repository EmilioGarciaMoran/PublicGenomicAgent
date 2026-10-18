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
    if not manifest.has_pedigree():
        raise ValueError(
            f"CaseManifest '{manifest.case_id}' no tiene pedigrí completo "
            "(pedigree + proband). No se puede aplicar análisis mendeliano."
        )
    return list(manifest.pedigree)


def case_manifest_to_rois(manifest: CaseManifest) -> list[GenomicRange]:
    if not manifest.has_rois():
        raise ValueError(
            f"CaseManifest '{manifest.case_id}' no tiene candidate_rois. "
            "El sistema debe derivarlos (por HPO o por búsqueda) antes "
            "de invocar tools que los requieran."
        )
    return list(manifest.candidate_rois)


def case_manifest_to_hpo_terms(manifest: CaseManifest) -> dict[str, list[str]]:
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
    return {
        "fetch_roi": manifest.has_rois(),
        "call_variants": manifest.has_rois(),
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

# Strings que el LLM a veces emite como valores en lugar de null real
_NULL_LIKE = {"null", "none", "", "n/a", "na", "-"}


def _is_null_like(value: object) -> bool:
    """True si el valor es un string que representa 'null' o similar."""
    if value is None:
        return True
    if isinstance(value, str) and value.strip().lower() in _NULL_LIKE:
        return True
    return False


def _extract_json_block(raw: str) -> dict:
    """Extrae el primer bloque JSON válido del output del LLM."""
    text = raw.strip()
    if text.startswith("```"):
        lines = text.splitlines()
        text = "\n".join(lines[1:])
        if text.rstrip().endswith("```"):
            text = text.rstrip()[:-3]

    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

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
    """Verifica que los samples del manifiesto aparecen en el texto."""
    text_lower = text.lower()
    verified: list[str] = []
    unverified: list[str] = []


    SYNONYMS = {
        "father": [
            "father", "padre", "papá", "the father", "el padre",
            "parents", "padres", "progenitors", "progenitores",
        ],
        "mother": [
            "mother", "madre", "mamá", "the mother", "la madre",
            "parents", "padres", "progenitors", "progenitoras",
        ],
        "proband": [
            "proband", "paciente", "patient", "hijo", "hija",
            "child", "son", "daughter", "male patient", "female patient",
            "the patient", "el paciente", "la paciente", "affected child",
        ],
        "son": ["son", "hijo", "male child"],
        "daughter": ["daughter", "hija", "female child"],
        "sibling": [
            "sibling", "hermano", "hermana", "brother", "sister",
            "older brother", "older sister",
        ],
        "brother": ["brother", "hermano", "older brother"],
        "sister": ["sister", "hermana", "older sister"],
    }
    for ind in manifest.pedigree:
        sample = ind.sample
        if sample.lower() in text_lower:
            verified.append(sample)
            continue
        synonyms = SYNONYMS.get(sample.lower(), [])
        if any(syn in text_lower for syn in synonyms):
            verified.append(sample)
            continue
        unverified.append(sample)

    return verified, unverified


def _sanitize_pedigree(
    manifest: CaseManifest,
) -> tuple[CaseManifest, list[str]]:
    """Limpia errores comunes del LLM en el pedigrí.

    1. Elimina samples con nombre "null", "none", "" o similar.
    2. Convierte referencias "father": "null" (string) a None.
    3. Elimina autorreferencias: un sample no puede ser su propio
       padre/madre.
    """
    warnings: list[str] = []
    cleaned: list[IndividualSpec] = []

    for ind in manifest.pedigree:
        # 1. Sample inválido
        if _is_null_like(ind.sample):
            warnings.append(
                f"Eliminado sample con nombre inválido: {ind.sample!r}."
            )
            continue

        # 2. Referencias father/mother con "null" string
        father = None if _is_null_like(ind.father) else ind.father
        mother = None if _is_null_like(ind.mother) else ind.mother

        # 3. Autorreferencias
        if father == ind.sample:
            warnings.append(
                f"Sample '{ind.sample}' se referenciaba a sí mismo como "
                f"padre; limpiado a null."
            )
            father = None
        if mother == ind.sample:
            warnings.append(
                f"Sample '{ind.sample}' se referenciaba a sí mismo como "
                f"madre; limpiado a null."
            )
            mother = None

        cleaned.append(
            IndividualSpec(
                sample=ind.sample,
                sex=ind.sex,
                affected=ind.affected,
                father=father,
                mother=mother,
            )
        )

    manifest.pedigree = cleaned
    return manifest, warnings


def _clean_affected_samples(
    manifest: CaseManifest,
) -> tuple[CaseManifest, list[str]]:
    """Recalcula `affected_samples` desde el pedigree.

    Fuente única de verdad: un sample está afectado si y solo si
    tiene `affected=True` en el pedigree. Esto elimina cualquier
    discrepancia introducida por el LLM.
    """
    warnings: list[str] = []
    original = list(manifest.affected_samples)

    derived = [ind.sample for ind in manifest.pedigree if ind.affected]

    # Detectar discrepancias para el warning
    if set(original) != set(derived):
        warnings.append(
            f"affected_samples recalculado desde pedigree: "
            f"LLM dijo {original}, pedigree dice {derived}."
        )

    manifest.affected_samples = derived
    return manifest, warnings


def _enforce_pedigree_consistency(
    manifest: CaseManifest,
) -> tuple[CaseManifest, list[str]]:
    """Añade samples referenciados como father/mother que no están
    definidos en el pedigree."""
    warnings: list[str] = []
    defined = {ind.sample for ind in manifest.pedigree}

    to_add: set[str] = set()
    for ind in manifest.pedigree:
        if ind.father and ind.father not in defined:
            to_add.add(ind.father)
        if ind.mother and ind.mother not in defined:
            to_add.add(ind.mother)

    for sample in sorted(to_add):
        manifest.pedigree.append(
            IndividualSpec(
                sample=sample,
                sex="U",
                affected=False,
                father=None,
                mother=None,
            )
        )
        warnings.append(
            f"Sample '{sample}' referenciado como padre/madre pero no "
            f"definido en el pedigrí; añadido con sex=U, affected=false."
        )

    return manifest, warnings
def _normalize_individual_fields(
    data: dict,
) -> tuple[dict, list[str]]:
    """Normaliza campos problemáticos en el pedigrí antes de Pydantic.

    El LLM a veces emite `"affected": null` en lugar de `true`/`false`.
    Pydantic lo rechaza. Convertimos:
      - affected: null -> false
      - sex: null o inválido -> "U"
    """
    warnings: list[str] = []

    pedigree = data.get("pedigree", [])
    if not isinstance(pedigree, list):
        return data, warnings

    for i, ind in enumerate(pedigree):
        if not isinstance(ind, dict):
            continue

        # affected: null o no-boolean -> false
        if ind.get("affected") is None:
            ind["affected"] = False
            warnings.append(
                f"pedigree[{i}] ({ind.get('sample', '?')}): "
                f"affected=null -> false."
            )
        elif not isinstance(ind.get("affected"), bool):
            ind["affected"] = bool(ind.get("affected"))
            warnings.append(
                f"pedigree[{i}] ({ind.get('sample', '?')}): "
                f"affected no booleano -> {ind['affected']}."
            )

        # sex: default "U"
        sex = ind.get("sex")
        if sex not in ("M", "F", "U"):
            ind["sex"] = "U"
            if sex not in (None, "U"):
                warnings.append(
                    f"pedigree[{i}] ({ind.get('sample', '?')}): "
                    f"sex inválido {sex!r} -> 'U'."
                )

    return data, warnings

def _sanitize_lists(
    data: dict,
) -> tuple[dict, list[str]]:
    """Corta listas sospechosamente largas antes de Pydantic.

    El LLM a veces alucina listas gigantes de HPOs, genes o ROIs
    inventados. Si una lista supera MAX_LIST_LEN, la vaciamos.
    """
    MAX_LIST_LEN = 10
    warnings: list[str] = []

    # hpo_terms: dict sample -> list[HPO]
    hpo_terms = data.get("hpo_terms")
    if isinstance(hpo_terms, dict):
        for sample, terms in list(hpo_terms.items()):
            if isinstance(terms, list) and len(terms) > MAX_LIST_LEN:
                hpo_terms[sample] = []
                warnings.append(
                    f"hpo_terms['{sample}'] tenía {len(terms)} términos "
                    f"(>{MAX_LIST_LEN}); probablemente alucinado, vaciado."
                )

    # candidate_rois y candidate_genes
    for key in ("candidate_rois", "candidate_genes"):
        val = data.get(key)
        if isinstance(val, list) and len(val) > MAX_LIST_LEN:
            data[key] = []
            warnings.append(
                f"{key} tenía {len(val)} elementos (>{MAX_LIST_LEN}); "
                f"probablemente alucinado, vaciado."
            )

    return data, warnings



def _validate_and_clean(
    data: dict,
    text: str,
) -> tuple[CaseManifest, list[str]]:
    """Valida y limpia un dict como CaseManifest.

    Pipeline:
      1. Pydantic valida el esquema (Nivel 1).
      2. _sanitize_pedigree: limpia "null", autorreferencias, etc.
      3. _clean_affected_samples: recalcula desde pedigree.
      4. _enforce_pedigree_consistency: completa referencias.
      5. _verify_samples_in_text: Nivel 2 contra el texto.
    """
    warnings: list[str] = []

    # Asegurar case_id
    if not data.get("case_id"):
        import uuid
        data["case_id"] = f"nl_case_{uuid.uuid4().hex[:8]}"
        warnings.append(f"case_id derivado automáticamente: {data['case_id']}")

    # Asegurar extractor
    if not data.get("extractor"):
        data["extractor"] = "llm"

    # Limpiar campos null
    for k in ("pedigree", "hpo_terms", "phenotype_text",
              "candidate_rois", "candidate_genes", "affected_samples"):
        if data.get(k) is None:
            data[k] = [] if k != "hpo_terms" and k != "phenotype_text" else {}
    # Normalizar campos problemáticos del pedigrí antes de Pydantic
        # Normalizar campos problemáticos del pedigrí antes de Pydantic
    data, wn = _normalize_individual_fields(data)
    warnings.extend(wn)

    # Cortar listas sospechosamente largas (alucinaciones de HPOs, genes)
    data, wl = _sanitize_lists(data)
    warnings.extend(wl)

    # Nivel 1: Pydantic
    manifest = CaseManifest.model_validate(data)


    # Post-procesado en orden
    manifest, w1 = _sanitize_pedigree(manifest)
    warnings.extend(w1)

    manifest, w2 = _clean_affected_samples(manifest)
    warnings.extend(w2)

    manifest, w3 = _enforce_pedigree_consistency(manifest)
    warnings.extend(w3)

    # Nivel 2: verificación contra el texto
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
    """Extrae un CaseManifest de un relato clínico con un LLM local."""
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

    cfg = load_config()
    client = build_llm_client(cfg.llm)

    user = render_extract_manifest_user(text, case_id=inp.case_id)
    raw = client.complete(EXTRACT_MANIFEST_SYSTEM_PROMPT, user)

    data = _extract_json_block(raw)
    manifest, warnings = _validate_and_clean(data, text)

    output_dir = Path(inp.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / f"{manifest.case_id}.manifest.json"
       # Sobreescribir el extractor con el modelo real (el LLM a veces inventa)
    real_extractor = f"llm:{cfg.llm.model}"
    if manifest.extractor != real_extractor:
        manifest.extractor = real_extractor

    # Escribir a disco
    output_dir = Path(inp.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / f"{manifest.case_id}.manifest.json"
    save_case_manifest(manifest, manifest_path)
   

    verified, unverified = _verify_samples_in_text(manifest, text)

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
    # Sobreescribir el extractor con el modelo real (el LLM a veces inventa)
    real_extractor = f"llm:{cfg.llm.model}"
    if manifest.extractor != real_extractor:
        manifest.extractor = real_extractor
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
