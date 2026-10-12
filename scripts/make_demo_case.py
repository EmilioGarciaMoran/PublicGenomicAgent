#!/usr/bin/env python3
"""Genera el caso demo (demo_trio) de forma reproducible y determinista.

Dos usos:

  1. Como script (CI, demos):
       python scripts/make_demo_case.py
     Escribe referencia, BAMs y un `manifest.json` de inspección.

  2. Como módulo (tests, CI in-memory):
       from scripts.make_demo_case import build_demo_case
       manifest = build_demo_case()   # CaseManifest, sin tocar disco

El manifest es un dump derivado; la fuente de verdad es
`build_demo_case()`. Los BAMs llevan @RG con SM:<sample>, necesario
para variant calling multi-muestra. La seed es fija, así que dos
ejecuciones producen los mismos ficheros.

Uso CLI:
    python scripts/make_demo_case.py
    python scripts/make_demo_case.py --check-only
    python scripts/make_demo_case.py --dump-manifest
    python scripts/make_demo_case.py --seed 123
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import pysam

from publicgenomicagent.tools.base import (
    CaseManifest,
    GenomicRange,
    IndividualSpec,
)


# ---------------------------------------------------------------------------
# Constantes del caso
# ---------------------------------------------------------------------------

CASE_DIR = Path("cases/demo_trio")
REF_DIR = CASE_DIR / "ref"
BAMS_DIR = CASE_DIR / "bams"
REF_FA = REF_DIR / "test_ref.fa"
REF_FAI = REF_DIR / "test_ref.fa.fai"
MANIFEST = CASE_DIR / "manifest.json"

CASE_ID = "demo_trio"
CHROM = "chr2"
CHROM_LEN = 200_000
ROI_START = 40_000
ROI_END = 60_000
ROI_LABEL = "NPHP1_demo_region"
CANDIDATE_GENES = ["NPHP1"]

# (sample, sex, affected, father, mother)
SAMPLES: list[tuple[str, str, bool, str | None, str | None]] = [
    ("father", "M", False, None, None),
    ("mother", "F", False, None, None),
    ("proband", "U", True, "father", "mother"),
]
PROBAND = "proband"

READ_LEN = 150
COVERAGE = 30
DEFAULT_SEED = 42


# ---------------------------------------------------------------------------
# Manifest in-memory (fuente de verdad)
# ---------------------------------------------------------------------------

def build_demo_case() -> CaseManifest:
    """Devuelve el CaseManifest del caso demo sin tocar disco."""
    return CaseManifest(
        case_id=CASE_ID,
        source="clinical",
        pedigree=[
            IndividualSpec(
                sample=s, sex=sex, affected=aff, father=father, mother=mother
            )
            for (s, sex, aff, father, mother) in SAMPLES
        ],
        proband=PROBAND,
        consanguinity=False,
        affected_samples=[s for (s, _, aff, _, _) in SAMPLES if aff],
        candidate_rois=[
            GenomicRange(
                chrom=CHROM,
                start=ROI_START,
                end=ROI_END,
                label=ROI_LABEL,
            )
        ],
        candidate_genes=list(CANDIDATE_GENES),
        confidence=1.0,
        extractor="manual",
        confirmed_by=None,
    )


# ---------------------------------------------------------------------------
# Referencia sintética
# ---------------------------------------------------------------------------

def make_reference(path: Path) -> None:
    """FASTA sintético de CHROM_LEN bases (ACGT repetido) + .fai."""
    path.parent.mkdir(parents=True, exist_ok=True)
    line_len = 80
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f">{CHROM}\n")
        written = 0
        base = "ACGT"
        while written < CHROM_LEN:
            take = min(line_len, CHROM_LEN - written)
            line = "".join(base[(written + i) % 4] for i in range(take))
            fh.write(line + "\n")
            written += take

    # .fai: nombre, longitud, offset del primer base, longitud de línea, ancho
    with path.open("rb") as fh:
        header = fh.readline()
        offset = len(header)
    fai = path.with_suffix(path.suffix + ".fai")
    with fai.open("w", encoding="utf-8") as fh:
        fh.write(f"{CHROM}\t{CHROM_LEN}\t{offset}\t{line_len}\t{line_len + 1}\n")


# ---------------------------------------------------------------------------
# BAMs
# ---------------------------------------------------------------------------

def _random_read(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def make_bam(sample: str, out_path: Path, seed: int) -> None:
    """BAM con reads de READ_LEN pb en [ROI_START, ROI_END], con @RG SM:<sample>."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
        "RG": [
            {
                "ID": sample,
                "SM": sample,
                "PL": "ILLUMINA",
                "LB": f"lib_{sample}",
            }
        ],
    }

    region_len = ROI_END - ROI_START
    n_reads = max(1, (region_len * COVERAGE) // READ_LEN)

    # Semilla por sample: determinista y estable entre ejecuciones.
    sample_seed = (seed * 1_000_003) ^ (sum(ord(c) for c in sample))
    rng = random.Random(sample_seed)

    tmp = out_path.with_suffix(".tmp.bam")
    with pysam.AlignmentFile(str(tmp), "wb", header=header) as fh:
        for i in range(n_reads):
            offset = rng.randrange(0, max(1, region_len - READ_LEN))
            ref_start = ROI_START + offset
            seq = _random_read(rng, READ_LEN)

            seg = pysam.AlignedSegment()
            seg.query_name = f"sim_{sample}_{i:07d}"
            seg.query_sequence = seq
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = ref_start
            seg.mapping_quality = 60
            seg.cigar = ((0, READ_LEN),)
            seg.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            seg.set_tag("RG", sample)
            fh.write(seg)

    pysam.sort("-o", str(out_path), str(tmp))
    pysam.index(str(out_path))
    if tmp.exists():
        tmp.unlink()


# ---------------------------------------------------------------------------
# Verificación
# ---------------------------------------------------------------------------

def check() -> list[str]:
    """Comprueba que todo existe y que el manifest valida."""
    issues: list[str] = []

    if not MANIFEST.exists():
        issues.append(f"falta {MANIFEST}")
    if not REF_FA.exists():
        issues.append(f"falta {REF_FA}")
    if not REF_FAI.exists():
        issues.append(f"falta {REF_FAI}")

    for sample, *_ in SAMPLES:
        bam = BAMS_DIR / f"{sample}.bam"
        bai = BAMS_DIR / f"{sample}.bam.bai"
        if not bam.exists():
            issues.append(f"falta {bam}")
        if not bai.exists():
            issues.append(f"falta {bai}")

    if MANIFEST.exists():
        try:
            raw = json.loads(MANIFEST.read_text(encoding="utf-8"))
            m = CaseManifest.model_validate(raw)
            m.validate_minimum()
            print(
                f"manifest OK: case_id={m.case_id} "
                f"samples={len(m.pedigree)} rois={len(m.candidate_rois)}"
            )
        except Exception as e:  # noqa: BLE001
            issues.append(f"manifest inválido: {type(e).__name__}: {e}")

    return issues


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def write_all(seed: int) -> None:
    """Regenera referencia, BAMs y manifest en disco."""
    print("Generando referencia sintética...")
    make_reference(REF_FA)

    for sample, *_ in SAMPLES:
        print(f"Generando BAM de {sample}...")
        make_bam(sample, BAMS_DIR / f"{sample}.bam", seed=seed)

    print("Escribiendo manifest (dump de inspección)...")
    manifest = build_demo_case()
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(
        json.dumps(manifest.model_dump(mode="json"), indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="No regenera; solo verifica que los ficheros existen y son válidos.",
    )
    parser.add_argument(
        "--dump-manifest",
        action="store_true",
        help="Regenera solo cases/demo_trio/manifest.json desde build_demo_case().",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Semilla para la generación de reads (default: {DEFAULT_SEED}).",
    )
    args = parser.parse_args()

    if args.dump_manifest and not args.check_only:
        print("Escribiendo manifest (dump de inspección)...")
        manifest = build_demo_case()
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        MANIFEST.write_text(
            json.dumps(manifest.model_dump(mode="json"), indent=2) + "\n",
            encoding="utf-8",
        )
    elif not args.check_only:
        write_all(seed=args.seed)

    issues = check()
    if issues:
        print(f"✗ {len(issues)} problemas:")
        for i in issues:
            print(f"  · {i}")
        return 1

    print(f"✓ Caso demo listo en {CASE_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

