#!/usr/bin/env python3
"""Genera un BAM sintético con UNA variante aplicada (SNP, del o ins).

Uso:
    make_bam_with_variant.py <ref.fa> <chrom> <pos_1based> <ref_allele> <alt_allele> \
        <sample> <out.bam> [coverage] [genotype]

Ejemplo:
    make_bam_with_variant.py ref.fa chr2 153 TT T KSA001_MAT out.bam 30 1/1
"""
from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

SAMTOOLS = str(Path("~/.pga/envs/pga-hts/bin/samtools").expanduser())
READ_LENGTH = 150
ERROR_RATE = 0.002
SEED = 42


def read_fasta(path: Path) -> tuple[str, str]:
    contig = None
    parts: list[str] = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if contig is None:
                    contig = line[1:].split()[0]
                continue
            parts.append(line)
    if contig is None:
        raise ValueError(f"FASTA sin cabecera: {path}")
    return contig, "".join(parts).upper()


def introduce_errors(seq: str, rng: random.Random) -> tuple[str, str]:
    out_bases = []
    out_quals = []
    for base in seq:
        q = int(rng.gauss(35, 5))
        q = max(2, min(40, q))
        if rng.random() < ERROR_RATE:
            new_base = rng.choice([b for b in "ACGT" if b != base])
            out_bases.append(new_base)
            out_quals.append(chr(min(q, 20) + 33))
        else:
            out_bases.append(base)
            out_quals.append(chr(q + 33))
    return "".join(out_bases), "".join(out_quals)


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt):
        return "snp"
    if len(ref) > len(alt):
        return "del"
    return "ins"


def build_read(
    ref_seq: str,
    start_1based: int,
    read_len: int,
    var_pos_1based: int,
    var_ref: str,
    var_alt: str,
    carries_variant: bool,
    rng: random.Random,
) -> tuple[str, str] | None:
    """Devuelve (seq, cigar) o None si la lectura no es válida.

    start_1based: posición 1-based del primer base de la lectura.
    var_pos_1based: posición 1-based del primer base del REF de la variante.
    """
    ref_idx = var_pos_1based - 1  # 0-based
    start_idx = start_1based - 1  # 0-based
    var_len_ref = len(var_ref)
    var_len_alt = len(var_alt)

    # Rango que cubre la lectura en la referencia (excluyendo la parte delecionada)
    read_end_ref = start_idx + read_len  # posición exclusiva en ref
    if start_idx < 0 or read_end_ref > len(ref_seq):
        return None

    # ¿La lectura cubre la variante?
    covers_variant = (start_idx <= ref_idx) and (ref_idx + var_len_ref <= read_end_ref)

    if not covers_variant or not carries_variant:
        # Lectura normal, 100% match
        seq = ref_seq[start_idx:read_end_ref]
        cigar = f"{len(seq)}M"
        return seq, cigar

    # La lectura cubre la variante y porta el alelo alt
    vtype = variant_type(var_ref, var_alt)

    if vtype == "snp":
        seq = list(ref_seq[start_idx:read_end_ref])
        rel = ref_idx - start_idx
        seq[rel:rel + var_len_ref] = list(var_alt)
        return "".join(seq), f"{len(seq)}M"

    if vtype == "del":
        # REF más largo que ALT: la lectura omite las bases extra
        rel_start = ref_idx - start_idx
        # Segmento antes de la deleción
        before = ref_seq[start_idx:ref_idx + len(var_alt)]  # mantiene la base ancla (parte del ALT)
        # Segmento después de la deleción
        del_len = var_len_ref - var_len_alt
        after_start = ref_idx + var_len_ref
        after = ref_seq[after_start:read_end_ref]
        seq = before + after
        # CIGAR: [antes]M [del_len]D [después]M
        cigar = f"{len(before)}M{del_len}D{len(after)}M"
        return seq, cigar

    if vtype == "ins":
        # REF más corto que ALT: la lectura añade bases extra
        rel_start = ref_idx - start_idx
        ins_len = var_len_alt - var_len_ref
        # Antes: hasta el final del REF
        before = ref_seq[start_idx:ref_idx + var_len_ref]
        # Insertamos las bases extra del ALT que van más allá del REF
        inserted = var_alt[var_len_ref:]
        # Después: desde el final del REF
        after_start = ref_idx + var_len_ref
        after = ref_seq[after_start:read_end_ref]
        seq = before + inserted + after
        # CIGAR: [antes]M [ins_len]I [después]M
        cigar = f"{len(before)}M{ins_len}I{len(after)}M"
        return seq, cigar

    return None


def main() -> int:
    if len(sys.argv) < 8:
        print(__doc__, file=sys.stderr)
        return 2

    ref_fa = Path(sys.argv[1]).expanduser().resolve()
    chrom = sys.argv[2]
    pos = int(sys.argv[3])
    ref_allele = sys.argv[4].upper()
    alt_allele = sys.argv[5].upper()
    sample = sys.argv[6]
    out_bam = Path(sys.argv[7]).expanduser().resolve()
    coverage = int(sys.argv[8]) if len(sys.argv) > 8 else 30
    genotype = sys.argv[9] if len(sys.argv) > 9 else "0/1"

    contig, ref_seq = read_fasta(ref_fa)
    seq_len = len(ref_seq)
    ref_idx = pos - 1

    # Verificar REF
    actual_ref = ref_seq[ref_idx:ref_idx + len(ref_allele)]
    if actual_ref != ref_allele:
        print(f"ERROR: REF no coincide en {chrom}:{pos}", file=sys.stderr)
        print(f"  esperado: {ref_allele}", file=sys.stderr)
        print(f"  encontrado: {actual_ref}", file=sys.stderr)
        return 1

    vtype = variant_type(ref_allele, alt_allele)
    print(f"Generando BAM con {vtype}")
    print(f"  Referencia: {ref_fa.name} ({seq_len} bp)")
    print(f"  Variante: {chrom}:{pos} {ref_allele}>{alt_allele}")
    print(f"  Genotipo: {genotype}")
    print(f"  Sample: {sample}")

    # Interpretar genotipo
    gt = genotype.replace("|", "/").split("/")
    n_alt = sum(1 for a in gt if a == "1")
    prob_alt = n_alt / len(gt)  # 0.0, 0.5 o 1.0

    rng = random.Random(SEED)
    n_reads = (seq_len * coverage) // READ_LENGTH
    reads = []
    n_with_variant = 0

    for i in range(n_reads):
        carries = rng.random() < prob_alt
        max_start = max(1, seq_len - READ_LENGTH + 1)
        start = rng.randint(1, max_start)

        result = build_read(
            ref_seq, start, READ_LENGTH,
            pos, ref_allele, alt_allele,
            carries, rng,
        )
        if result is None:
            continue

        seq, cigar = result
        if carries and cigar != f"{len(seq)}M":
            n_with_variant += 1
        reads.append((f"r{i+1:08d}", start, seq, cigar))

    # Escribir SAM
    sam_path = out_bam.with_suffix(".sam")
    with sam_path.open("w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        f.write(f"@SQ\tSN:{contig}\tLN:{seq_len}\n")
        f.write(f"@RG\tID:rg1\tSM:{sample}\tPL:ILLUMINA\n")
        for name, start, seq, cigar in sorted(reads, key=lambda x: x[1]):
            noisy, quals = introduce_errors(seq, rng)
            f.write(
                f"{name}\t0\t{contig}\t{start}\t60\t{cigar}\t*\t0\t0\t"
                f"{noisy}\t{quals}\tRG:Z:rg1\n"
            )

    subprocess.run([SAMTOOLS, "view", "-b", "-o", str(out_bam), str(sam_path)], check=True)
    subprocess.run([SAMTOOLS, "index", str(out_bam)], check=True)
    sam_path.unlink()

    print(f"OK: {out_bam}")
    print(f"  Lecturas totales: {len(reads)}")
    print(f"  Lecturas con variante (CIGAR != full match): {n_with_variant}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
