"""Tests del generador del caso demo (siembra de variante).

Verifican que:
  - father  es heterocigoto (0/1) en VARIANT_POS
  - mother  es heterocigoto (0/1) en VARIANT_POS
  - proband es homocigoto alt (1/1) en VARIANT_POS
  - la base del proband difiere de la referencia

El test regenera el caso con `python scripts/make_demo_case.py` y luego
lee el pileup directamente con pysam. Es determinista (seed fija).
"""
from __future__ import annotations

import subprocess
from collections import Counter
from pathlib import Path

import pysam
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "make_demo_case.py"
CASE_DIR = REPO_ROOT / "cases" / "demo_trio"
BAMS_DIR = CASE_DIR / "bams"

# Debe coincidir con scripts/make_demo_case.py
CHROM = "chr2"
VARIANT_POS_1BASED = 50_001
VARIANT_POS_0BASED = VARIANT_POS_1BASED - 1


@pytest.fixture(scope="module")
def demo_case() -> Path:
    """Regenera el caso demo antes de ejecutar los tests del módulo."""
    result = subprocess.run(
        ["python", str(SCRIPT)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"make_demo_case.py falló:\nSTDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    return CASE_DIR


def _pileup_at(bam_path: Path, chrom: str, pos_0based: int) -> dict[str, int]:
    """Devuelve {base: count} para la posición 0-based indicada."""
    counts: Counter[str] = Counter()
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for col in bam.pileup(
            chrom, pos_0based, pos_0based + 1, truncate=True, stepper="nofilter"
        ):
            for p in col.pileups:
                if p.query_position is not None:
                    base = p.alignment.query_sequence[p.query_position]
                    counts[base] += 1
    return dict(counts)


def _non_zero(counts: dict[str, int]) -> dict[str, int]:
    return {b: c for b, c in counts.items() if c > 0}


def test_father_is_heterozygous(demo_case: Path) -> None:
    counts = _pileup_at(BAMS_DIR / "father.bam", CHROM, VARIANT_POS_0BASED)
    nz = _non_zero(counts)
    assert len(nz) == 2, (
        f"father debería ser HET en {VARIANT_POS_1BASED}, pileup={counts}"
    )


def test_mother_is_heterozygous(demo_case: Path) -> None:
    counts = _pileup_at(BAMS_DIR / "mother.bam", CHROM, VARIANT_POS_0BASED)
    nz = _non_zero(counts)
    assert len(nz) == 2, (
        f"mother debería ser HET en {VARIANT_POS_1BASED}, pileup={counts}"
    )


def test_proband_is_homozygous_alt(demo_case: Path) -> None:
    counts = _pileup_at(BAMS_DIR / "proband.bam", CHROM, VARIANT_POS_0BASED)
    nz = _non_zero(counts)
    assert len(nz) == 1, (
        f"proband debería ser HOM ALT en {VARIANT_POS_1BASED}, pileup={counts}"
    )


def test_proband_base_differs_from_reference(demo_case: Path) -> None:
    """La base del proband debe ser distinta de la referencia en esa posición."""
    ref_path = CASE_DIR / "ref" / "test_ref.fa"
    with pysam.FastaFile(str(ref_path)) as fasta:
        ref_base = fasta.fetch(CHROM, VARIANT_POS_0BASED, VARIANT_POS_0BASED + 1)

    counts = _pileup_at(BAMS_DIR / "proband.bam", CHROM, VARIANT_POS_0BASED)
    nz = _non_zero(counts)
    assert len(nz) == 1, f"proband debería ser HOM, pileup={counts}"
    proband_base = next(iter(nz))

    assert proband_base != ref_base, (
        f"proband {proband_base} == referencia {ref_base}; "
        f"la variante no se ha sembrado"
    )

