#!/usr/bin/env python3
"""Genera el caso demo (demo_trio) de forma reproducible y determinista.

Produce:
  cases/demo_trio/ref/test_ref.fa       (chr2 sintético, 200 kb)
  cases/demo_trio/ref/test_ref.fa.fai
  cases/demo_trio/bams/{father,mother,proband}.bam(+.bai)
  cases/demo_trio/manifest.json

Idempotente y sin dependencias externas: solo pysam.
Los BAMs llevan @RG con SM:<sample>, necesario para variant calling
multi-muestra. La seed es fija, así que dos ejecuciones producen
exactamente los mismos ficheros.

Los reads se generan a partir de la referencia real (con una única
variante sembrada en VARIANT_POS, 1-based) más una tasa pequeña de
errores de secuenciación. La zygosidad del sample determina la mezcla
de reads del haplotipo de referencia y del haplotipo mutado:
  - father:  0/1 (50% ref, 50% alt)
  - mother:  0/1 (50% ref, 50% alt)
  - proband: 1/1 (100% alt)
Así mendelian_filter detecta la variante como auto_rec_hom.

Uso:
    python scripts/make_demo_case.py
    python scripts/make_demo_case.py --check-only   # solo verifica
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
# Siembra de variante
# ---------------------------------------------------------------------------

# Variante sembrada: SNV en VARIANT_POS (1-based, dentro de la ROI).
# La base de referencia se deriva del FASTA en tiempo de ejecución;
# VARIANT_ALT se elige distinta para garantizar una sustitución real.
VARIANT_POS = 50_001
VARIANT_ALT = "C"

# Zygosidad por sample: probabilidad de que un read venga del haplotipo alt.
ZYGOSITY: dict[str, float] = {
    "father": 0.5,   # 0/1
    "mother": 0.5,   # 0/1
    "proband": 1.0,  # 1/1
}

# Tasa de error de secuenciación por base (sustituciones). Da realismo sin
# ahogar la señal de la variante sembrada.
ERROR_RATE = 0.001


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

def _load_reference(ref_path: Path) -> pysam.FastaFile:
    """Abre el FASTA de referencia; asegura el .fai antes."""
    if not ref_path.exists():
        raise FileNotFoundError(f"referencia no encontrada: {ref_path}")
    if not ref_path.with_suffix(ref_path.suffix + ".fai").exists():
        pysam.faidx(str(ref_path))
    return pysam.FastaFile(str(ref_path))


def _build_haplotypes(fasta: pysam.FastaFile) -> tuple[str, str, str]:
    """Construye (ref_haplotype, alt_haplotype, ref_base_en_VARIANT_POS).

    El haplotipo cubre [ROI_START, ROI_END] y en la posición VARIANT_POS
    sustituye la base de referencia por VARIANT_ALT.
    """
    ref_seq = fasta.fetch(CHROM, ROI_START, ROI_END).upper()
    # VARIANT_POS es 1-based; convertimos a 0-based dentro de la ROI.
    rel_pos = (VARIANT_POS - 1) - ROI_START
    if rel_pos < 0 or rel_pos >= len(ref_seq):
        raise ValueError(
            f"VARIANT_POS={VARIANT_POS} fuera de la ROI "
            f"[{ROI_START}, {ROI_END})"
        )
    ref_base = ref_seq[rel_pos]
    alt_seq = ref_seq[:rel_pos] + VARIANT_ALT + ref_seq[rel_pos + 1:]
    return ref_seq, alt_seq, ref_base


def _mutate_read(seq: str, rng: random.Random) -> tuple[str, str]:
    """Aplica errores aleatorios con probabilidad ERROR_RATE.

    Devuelve (seq_mutada, qual_string). Bases sin error: Q40 ('I').
    Bases con error: Q10 ('+'), para que el caller las penalice.
    """
    bases = list(seq)
    quals = ["I"] * len(seq)
    for i in range(len(bases)):
        if rng.random() < ERROR_RATE:
            orig = bases[i]
            choices = [b for b in "ACGT" if b != orig]
            bases[i] = rng.choice(choices)
            quals[i] = "+"
    return "".join(bases), "".join(quals)


def make_bam(sample: str, out_path: Path, seed: int, ref_path: Path) -> None:
    """BAM con reads de READ_LEN pb en [ROI_START, ROI_END], con @RG SM:<sample>.

    Los reads se generan a partir de la referencia real (con una única
    variante sembrada en VARIANT_POS) más una tasa pequeña de errores
    de secuenciación. La zygosidad del sample determina la mezcla de
    reads del haplotipo de referencia y del haplotipo mutado.
    """
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

    fasta = _load_reference(ref_path)
    try:
        ref_hap, alt_hap, ref_base = _build_haplotypes(fasta)
    finally:
        fasta.close()

    if ref_base == VARIANT_ALT:
        raise ValueError(
            f"VARIANT_REF coincide con VARIANT_ALT en {CHROM}:{VARIANT_POS} "
            f"({ref_base}); elige otra posición o ALT"
        )

    zygosity = ZYGOSITY.get(sample, 0.0)
    region_len = ROI_END - ROI_START
    n_reads = max(1, (region_len * COVERAGE) // READ_LEN)

    sample_seed = (seed * 1_000_003) ^ sum(ord(c) for c in sample)
    rng = random.Random(sample_seed)

    tmp = out_path.with_suffix(".tmp.bam")
    with pysam.AlignmentFile(str(tmp), "wb", header=header) as fh:
        for i in range(n_reads):
            use_alt = rng.random() < zygosity
            hap = alt_hap if use_alt else ref_hap

            offset = rng.randrange(0, max(1, region_len - READ_LEN))
            read_seq = hap[offset: offset + READ_LEN]

            read_seq, qual_str = _mutate_read(read_seq, rng)

            seg = pysam.AlignedSegment()
            seg.query_name = f"sim_{sample}_{i:07d}"
            seg.query_sequence = read_seq
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = ROI_START + offset
            seg.mapping_quality = 60
            seg.cigar = ((0, READ_LEN),)
            seg.query_qualities = pysam.qualitystring_to_array(qual_str)
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
        make_bam(sample, BAMS_DIR / f"{sample}.bam", seed=seed, ref_path=REF_FA)

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
