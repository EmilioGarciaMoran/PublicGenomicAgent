#!/usr/bin/env python3
"""Genera un BAM sintético de un paciente desde su VCF.

Enfoque: cobertura uniforme a lo largo del ROI. Cada lectura se
construye aplicando el haplotipo completo del paciente (todas sus
variantes) en el tramo que cubre.

Uso:
    make_patient_bam.py <vcf.gz> <reference.fa> <sample> <out.bam> [coverage]
"""
from __future__ import annotations

import gzip
import random
import subprocess
import sys
from pathlib import Path

try:
    import publicgenomicagent  # noqa: F401
except ImportError:
    raise SystemExit(
        "ERROR: este script requiere el entorno pga-core.\n"
        "Ejecuta con:\n"
        "  ~/.pga/envs/pga-core/bin/python " + __file__
    )


SAMTOOLS = str(Path("~/.pga/envs/pga-hts/bin/samtools").expanduser())
READ_LENGTH = 150
SEED = 1234

# Parámetros de simulación realista
ERROR_RATE = 0.002          # 0.2% por base
TRANSITION_BIAS = 0.7       # 70% de errores son transiciones (A<->G, C<->T)
QUAL_MEAN = 35              # calidad media Q35 (~0.03% error esperado)
QUAL_STD = 5
BASES = "ACGT"
TRANSITIONS = {"A": "G", "G": "A", "C": "T", "T": "C"}
TRANSVERSIONS = {
    "A": ["C", "T"],
    "C": ["A", "G"],
    "G": ["C", "T"],
    "T": ["A", "G"],
}


def read_reference(fasta: Path) -> tuple[str, str]:
    contig = None
    parts: list[str] = []
    with fasta.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if contig is None:
                    contig = line[1:].split()[0]
                continue
            parts.append(line)
    if contig is None:
        raise ValueError(f"FASTA sin cabecera: {fasta}")
    return contig, "".join(parts)


def read_vcf_sample(vcf_gz: Path, sample: str) -> list[dict]:
    opener = gzip.open if str(vcf_gz).endswith(".gz") else open
    variants: list[dict] = []
    samples: list[str] = []
    sample_idx = -1

    with opener(vcf_gz, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                if sample not in samples:
                    raise ValueError(
                        f"sample '{sample}' no está en {vcf_gz}. "
                        f"Disponibles: {samples}"
                    )
                sample_idx = samples.index(sample)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            if len(ref) != 1 or len(alt) != 1:
                continue

            fmt = fields[8].split(":")
            if "GT" not in fmt:
                continue
            gt_idx = fmt.index("GT")
            sample_field = fields[9 + sample_idx].split(":")
            if gt_idx >= len(sample_field):
                continue
            gt = sample_field[gt_idx]
            if gt in ("./.", "."):
                continue

            variants.append({
                "pos": pos,
                "ref": ref.upper(),
                "alt": alt.upper(),
                "gt": gt,
            })

    return variants


def genotype_to_alleles(gt: str) -> tuple[int, int]:
    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)
    return (int(parts[0]), int(parts[1]))


def build_haplotype_sequence(
    ref_seq: str,
    variants: list[dict],
    allele_idx: int,
) -> str:
    """Construye la secuencia haploide completa aplicando, para cada
    variante, el alelo correspondiente al haplotype `allele_idx` (0 o 1).
    """
    hap = list(ref_seq)
    for v in variants:
        alleles = genotype_to_alleles(v["gt"])
        allele = alleles[allele_idx]
        if allele == 1:
            idx = v["pos"] - 1
            hap[idx] = v["alt"]
    return "".join(hap)


def introduce_errors(seq: str, rng: random.Random) -> tuple[str, str]:
    """Introduce errores de secuenciación y devuelve (seq_con_errores, qual_string).

    - ERROR_RATE: probabilidad de error por base.
    - TRANSITION_BIAS: fracción de errores que son transiciones.
    - Calidad variable alrededor de QUAL_MEAN con distribución normal.
    """
    out_bases = []
    out_quals = []
    for base in seq:
        q = int(rng.gauss(QUAL_MEAN, QUAL_STD))
        q = max(2, min(40, q))
        if rng.random() < ERROR_RATE:
            # Error: transición o transversión
            if rng.random() < TRANSITION_BIAS and base in TRANSITIONS:
                new_base = TRANSITIONS[base]
            elif base in TRANSVERSIONS:
                new_base = rng.choice(TRANSVERSIONS[base])
            else:
                new_base = rng.choice([b for b in BASES if b != base])
            out_bases.append(new_base)
            # Calidad baja en la base errónea
            out_quals.append(chr(min(q, 20) + 33))
        else:
            out_bases.append(base)
            out_quals.append(chr(q + 33))
    return "".join(out_bases), "".join(out_quals)


def main():
    if len(sys.argv) < 5:
        print(__doc__, file=sys.stderr)
        sys.exit(2)

    vcf_gz = Path(sys.argv[1]).expanduser().resolve()
    ref_fa = Path(sys.argv[2]).expanduser().resolve()
    sample = sys.argv[3]
    out_bam = Path(sys.argv[4]).expanduser().resolve()
    coverage = int(sys.argv[5]) if len(sys.argv) > 5 else 30

    contig, ref_seq = read_reference(ref_fa)
    variants = read_vcf_sample(vcf_gz, sample)

    rng = random.Random(SEED)
    seq_len = len(ref_seq)

    # Construir los dos haplotipos del paciente
    hap0 = build_haplotype_sequence(ref_seq, variants, 0)
    hap1 = build_haplotype_sequence(ref_seq, variants, 1)

    # Generar lecturas con cobertura uniforme
    # Cada lectura empieza en una posición aleatoria y viene de un haplotipo aleatorio
    n_reads = (seq_len * coverage) // READ_LENGTH
    reads: list[tuple[str, int, str]] = []

    for i in range(n_reads):
        # Elegir haplotipo aleatoriamente (50/50)
        hap = hap0 if rng.random() < 0.5 else hap1

        # Posición de inicio aleatoria (1-based)
        start = rng.randint(1, max(1, seq_len - READ_LENGTH + 1))
        end = min(seq_len, start + READ_LENGTH - 1)
        seq = hap[start - 1:end]
        length = len(seq)

        name = f"r{i+1:08d}"
        reads.append((name, start, seq))

    # Escribir SAM ordenado por posición
    sam_path = out_bam.with_suffix(".sam")
    with sam_path.open("w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        f.write(f"@SQ\tSN:{contig}\tLN:{seq_len}\n")
        f.write(f"@RG\tID:rg1\tSM:{sample}\tPL:ILLUMINA\n")

        for name, start, seq in sorted(reads, key=lambda x: x[1]):
            noisy_seq, quals = introduce_errors(seq, rng)
            f.write(
                f"{name}\t0\t{contig}\t{start}\t60\t{len(seq)}M\t*\t0\t0\t"
                f"{noisy_seq}\t{quals}\tRG:Z:rg1\n"
            )

    subprocess.run(
        [SAMTOOLS, "view", "-b", "-o", str(out_bam), str(sam_path)],
        check=True,
    )
    subprocess.run([SAMTOOLS, "index", str(out_bam)], check=True)
    sam_path.unlink()

    print(f"OK: {out_bam}")
    print(f"Lecturas: {len(reads)}  Sample: {sample}  Contig: {contig}  Cobertura objetivo: {coverage}x")


if __name__ == "__main__":
    main()
