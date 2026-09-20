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
