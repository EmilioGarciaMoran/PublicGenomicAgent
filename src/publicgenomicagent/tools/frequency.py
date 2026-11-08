"""annotate_frequency: frecuencia poblacional gnomAD v4 (global y Oriente Medio).

Añade al VCF los campos INFO PGA_AF, PGA_AF_MID, PGA_AC, PGA_AN y
PGA_AN_MID para cada variante encontrada en gnomAD v4 (exoma + genoma).
Una variante ausente en gnomAD no recibe campos: `mendelian_filter` la
trata como no observada (rara). La cabecera declara fuente y fecha, de
modo que el filtro puede exigir que la anotación se hizo de verdad.

Consulta solo las regiones de interés (ROIs), en ventanas de <= 50 kb,
con caché en ~/.pga/cache/gnomad_af: tras la primera consulta el
análisis es reproducible y funciona sin red.
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

from ..knowledge.gnomad import fetch_region_af, variant_frequencies
from .base import AnnotateFrequencyInput, AnnotateFrequencyOutput
from .mendelian import _load_vcf, _write_subvcf

HEADER_MARK = "##PGA_frequency_source="
INFO_LINES = [
    '##INFO=<ID=PGA_AF,Number=A,Type=Float,Description="gnomAD v4 AF (exoma+genoma)">',
    '##INFO=<ID=PGA_AF_MID,Number=A,Type=Float,Description="gnomAD v4 AF Oriente Medio (mid)">',
    '##INFO=<ID=PGA_AC,Number=A,Type=Integer,Description="gnomAD v4 AC (exoma+genoma)">',
    '##INFO=<ID=PGA_AN,Number=1,Type=Integer,Description="gnomAD v4 AN (exoma+genoma)">',
    '##INFO=<ID=PGA_AN_MID,Number=1,Type=Integer,Description="gnomAD v4 AN Oriente Medio">',
    '##INFO=<ID=PGA_REF_MISMATCH,Number=0,Type=Flag,Description="REF distinto del de gnomAD '
    'en esa posición: la referencia usada no es GRCh38; frecuencia no aplicable">',
]


def annotate_lines(header: list[str], lines: list[str], lookup: dict) -> tuple[list[str], list[str], dict]:
    """Añade los INFO a cada línea según `lookup[(chrom, pos, ref, alt)] -> freqs`."""
    stats = {"n": 0, "found": 0, "absent": 0, "multiallelic_skipped": 0, "ref_mismatch": 0}
    refs_at: dict[tuple[str, int], set[str]] = {}
    for (c, p, r, _a) in lookup:
        refs_at.setdefault((c, p), set()).add(r)
    out = []
    for line in lines:
        f = line.rstrip("\n").split("\t")
        stats["n"] += 1
        alts = f[4].split(",")
        if len(alts) != 1:
            stats["multiallelic_skipped"] += 1
            out.append(line)
            continue
        chrom, pos = f[0].removeprefix("chr"), int(f[1])
        fr = lookup.get((chrom, pos, f[3], alts[0]))
        if fr is None and refs_at.get((chrom, pos)) and f[3] not in refs_at[(chrom, pos)]:
            stats["ref_mismatch"] += 1
            f[7] = ";".join(([] if f[7] in (".", "") else [f[7]]) + ["PGA_REF_MISMATCH"])
            out.append("\t".join(f) + "\n")
            continue
        if fr is None or fr["an"] == 0:
            stats["absent"] += 1
            out.append(line)
            continue
        stats["found"] += 1
        add = [f"PGA_AF={fr['af']:.6g}", f"PGA_AC={fr['ac']}", f"PGA_AN={fr['an']}"]
        if fr["af_mid"] is not None:
            add += [f"PGA_AF_MID={fr['af_mid']:.6g}", f"PGA_AN_MID={fr['an_mid']}"]
        f[7] = ";".join(([] if f[7] in (".", "") else [f[7]]) + add)
        out.append("\t".join(f) + "\n")
    new_header = [h for h in header if not h.startswith("#CHROM")]
    new_header += [l + "\n" for l in INFO_LINES]
    new_header.append(f"{HEADER_MARK}gnomad_r4;date={dt.date.today().isoformat()}\n")
    new_header += [h for h in header if h.startswith("#CHROM")]
    return new_header, out, stats


def annotate_frequency(inp: AnnotateFrequencyInput) -> AnnotateFrequencyOutput:
    header, data, _ = _load_vcf(Path(inp.vcf))
    lookup: dict = {}
    methods: set[str] = set()
    for r in inp.regions:
        for v in fetch_region_af(r.chrom, r.start, r.end, refresh=inp.refresh):
            fr = variant_frequencies(v)
            methods.add(fr["method"])
            lookup[(r.chrom.removeprefix("chr"), int(v["pos"]), v["ref"], v["alt"])] = fr
    new_header, lines, stats = annotate_lines(header, data, lookup)
    out = Path(inp.output_vcf)
    out.parent.mkdir(parents=True, exist_ok=True)
    _write_subvcf(out, new_header, lines)
    return AnnotateFrequencyOutput(
        ok=True, tool="annotate_frequency",
        message=(f"gnomAD v4: {stats['found']}/{stats['n']} variantes con frecuencia, "
                 f"{stats['absent']} ausentes"
                 + (f"; AVISO: {stats['ref_mismatch']} con REF distinto de gnomAD "
                    "(¿referencia que no es GRCh38?)" if stats["ref_mismatch"] else "")),
        output_vcf=out, n_variants=stats["n"], n_found=stats["found"],
        n_absent=stats["absent"], n_multiallelic_skipped=stats["multiallelic_skipped"],
        n_ref_mismatch=stats["ref_mismatch"],
        gnomad_variants_in_regions=len(lookup), method=",".join(sorted(methods)) or "sin_datos",
        outputs={"vcf": str(out)},
    )
