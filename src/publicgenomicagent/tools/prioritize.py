"""Variant prioritisation by simple clinical scoring rules.

This is NOT a full ACMG implementation. It is a pragmatic ranking
that combines:
  - ClinVar significance (CLNSIG)
  - Proband genotype (1/1 gains points)
  - QUAL threshold
  - INDEL flag

The output is a VCF with the same variants (unchanged) plus a
sorted ranking list. The idea is to give the clinician a
prioritised shortlist instead of a raw VCF.
"""
from __future__ import annotations

import gzip
from pathlib import Path

from .base import (
    PrioritizeVariantsInput,
    PrioritizeVariantsOutput,
    VariantRank,
)


def _score_variant(
    *,
    clnsig: str | None,
    gt_proband: str | None,
    qual: float | None,
    is_indel: bool,
) -> tuple[float, list[str]]:
    score = 0.0
    reasons: list[str] = []
    if clnsig:
        c = clnsig.lower()
        if "pathogenic" in c and "likely" not in c:
            score += 3
            reasons.append("ClinVar: Pathogenic")
        elif "likely_pathogenic" in c or "likely pathogenic" in c:
            score += 2
            reasons.append("ClinVar: Likely pathogenic")
        elif "benign" in c:
            score -= 3
            reasons.append("ClinVar: Benign")
        elif "uncertain" in c or "vus" in c:
            reasons.append("ClinVar: VUS")
    if gt_proband in ("1/1", "1|1"):
        score += 2
        reasons.append("proband hom alt")
    if qual is not None and qual > 100:
        score += 1
        reasons.append("QUAL>100")
    if is_indel:
        score += 1
        reasons.append("INDEL")
    return score, reasons


def prioritize_variants(
    inp: PrioritizeVariantsInput,
) -> PrioritizeVariantsOutput:
    """Lee un VCF, puntúa cada variante y devuelve un ranking.

    El VCF de salida es idéntico al de entrada (no reordena el
    fichero). El ranking va en el output Pydantic.
    """
    vcf = Path(inp.vcf).expanduser().resolve()
    out_vcf = Path(inp.output_vcf).expanduser().resolve()
    if not vcf.exists():
        raise FileNotFoundError(f"VCF no existe: {vcf}")

    # Copiar el VCF (nos limitamos a copiarlo, no reordenamos)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(vcf, "rb") as fin, gzip.open(out_vcf, "wb") as fout:
        fout.write(fin.read())

    # Parsear y puntuar
    samples: list[str] = []
    ranking: list[VariantRank] = []
    total = 0

    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            total += 1

            chrom, pos, _id, ref, alt, qual_s, _filt = fields[:7]
            info_raw = fields[7]
            info = {}
            for kv in info_raw.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info[k] = v

            # Genotipo del probando
            gt_proband = None
            if inp.proband and inp.proband in samples and len(fields) >= 10:
                fmt = fields[8].split(":")
                if "GT" in fmt:
                    gt_idx = fmt.index("GT")
                    s_idx = samples.index(inp.proband)
                    parts = fields[9 + s_idx].split(":")
                    if gt_idx < len(parts):
                        gt_proband = parts[gt_idx]

            # Genotipos de todos
            genotypes = {}
            if len(fields) >= 10 and samples:
                fmt = fields[8].split(":")
                gt_idx = fmt.index("GT") if "GT" in fmt else 0
                for i, s in enumerate(samples):
                    parts = fields[9 + i].split(":")
                    if gt_idx < len(parts):
                        genotypes[s] = parts[gt_idx]

            try:
                qual_val = float(qual_s) if qual_s != "." else None
            except ValueError:
                qual_val = None

            is_indel = info.get("INDEL") == "1" or "INDEL" in info

            score, reasons = _score_variant(
                clnsig=info.get("CLNSIG"),
                gt_proband=gt_proband,
                qual=qual_val,
                is_indel=is_indel,
            )

            ranking.append(VariantRank(
                rank=0,
                chrom=chrom,
                pos=int(pos),
                ref=ref,
                alt=alt,
                qual=qual_val,
                genotypes=genotypes,
                clnsig=info.get("CLNSIG"),
                clndn=info.get("CLNDN"),
                score=score,
                reasons=reasons,
            ))

    # Ordenar por score descendente, luego por QUAL
    ranking.sort(key=lambda v: (v.score, v.qual or 0), reverse=True)
    for i, v in enumerate(ranking):
        v.rank = i + 1

    top = ranking[: inp.top_n]

    return PrioritizeVariantsOutput(
        ok=True,
        tool="prioritize_variants",
        message=f"{total} variantes puntuadas, top {len(top)}",
        outputs={"vcf": str(out_vcf)},
        output_vcf=out_vcf,
        ranking=top,
        variants_total=total,
        variants_scored=len(ranking),
    )
