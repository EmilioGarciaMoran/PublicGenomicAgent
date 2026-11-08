#!/usr/bin/env python3
"""E2: caso realista con referencia GRCh38 real y una variante recesiva sembrada.

Construye, sin datos de pacientes:
  1. Referencia: el cromosoma con su longitud real de GRCh38, con la
     secuencia real de la ROI (UCSC) y N en el resto.
  2. Comprueba que el REF de cada SNV de gnomAD v4 en la ROI coincide con
     esa secuencia (valida referencia y coordenadas; aborta si no).
  3. Familia consanguínea (padres primos; 2 afectados, 1 sano; abuelos y
     tíos en el PED sin BAM). Los fundadores reciben alelos según la
     frecuencia real de gnomAD (Oriente Medio si existe, si no global), con
     recombinación al transmitir.
  4. Siembra una SNV rara real de gnomAD como recesiva: padres 0/1,
     afectados 1/1, sano 0/1.
  5. BAMs de los 5 secuenciados, PED y manifiesto con la variante sembrada.
  6. Ejecuta `pga run-family --max-af 0.01` y evalúa los criterios de E2:
     la sembrada está entre las candidatas recesivas, las comunes reales se
     filtran y las candidatas son pocas. Escribe e2_report.md.

Sin --relatedness-sites: una sola ROI no tiene sitios suficientes para la
verificación familiar (sería "no verificable").

Uso (con red la primera vez; después todo sale de caché):
  pga-env
  python scripts/make_e2_case.py --workdir ~/pga_e2
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

PED = """FAM GF 0 0 1 1
FAM GM 0 0 2 1
FAM A GF GM 1 1
FAM B GF GM 2 1
FAM SA 0 0 2 1
FAM SB 0 0 1 1
FAM FA A SA 1 1
FAM MO SB B 2 1
FAM P1 FA MO 1 2
FAM P2 FA MO 2 2
FAM S1 FA MO 1 1
"""
PARENTS = {"GF": (None, None), "GM": (None, None), "A": ("GF", "GM"), "B": ("GF", "GM"),
           "SA": (None, None), "SB": (None, None), "FA": ("A", "SA"), "MO": ("SB", "B"),
           "P1": ("FA", "MO"), "P2": ("FA", "MO"), "S1": ("FA", "MO")}
SEQ = ["FA", "MO", "P1", "P2", "S1"]
SEEDED_GT = {"FA": 1, "MO": 1, "P1": 2, "P2": 2, "S1": 1}
GRCH38_LEN = {"chr22": 50_818_468, "chr2": 242_193_529, "chr1": 248_956_422}


def log(m: str) -> None:
    print(f"[e2] {m}", flush=True)


def read_slice(path: Path) -> str:
    return "".join(l.strip() for l in path.read_text().splitlines() if not l.startswith(">")).upper()


def build_reference(work: Path, chrom: str, start: int, seq: str, length: int) -> Path:
    import pysam
    ref = work / "ref" / f"grch38_{chrom}_roi.fa"
    if ref.exists() and Path(str(ref) + ".fai").exists():
        return ref
    ref.parent.mkdir(parents=True, exist_ok=True)
    rest = length - (start - 1 + len(seq))
    segments = [("N", start - 1), (seq, None), ("N", rest)]
    col = 0
    with ref.open("w") as f:
        f.write(f">{chrom}\n")
        for content, count in segments:
            if count is not None:                       # bloque de N, por trozos
                remaining = count
                while remaining:
                    take = min(remaining, 60 - col)
                    f.write("N" * take); col += take; remaining -= take
                    if col == 60:
                        f.write("\n"); col = 0
                    if col == 0 and remaining >= 60 * 100_000:   # rápido: líneas completas
                        lines = remaining // 60
                        f.write(("N" * 60 + "\n") * lines); remaining -= lines * 60
            else:
                i = 0
                while i < len(content):
                    take = min(len(content) - i, 60 - col)
                    f.write(content[i:i + take]); col += take; i += take
                    if col == 60:
                        f.write("\n"); col = 0
        if col:
            f.write("\n")
    pysam.faidx(str(ref))
    return ref


def snvs_matching_reference(variants: list[dict], seq: str, start: int):
    from publicgenomicagent.knowledge.gnomad import variant_frequencies
    ok, mism = [], 0
    for v in variants:
        if len(v["ref"]) != 1 or len(v["alt"]) != 1:
            continue
        p = int(v["pos"])
        if not start <= p < start + len(seq):
            continue
        if seq[p - start] != v["ref"]:
            mism += 1
            continue
        fr = variant_frequencies(v)
        if fr["an"]:
            ok.append((p, v["ref"], v["alt"], fr))
    return ok, mism


def simulate(sites, seed_idx: int, rng: random.Random, recomb: float = 0.02):
    def af(fr):
        return fr["af_mid"] if fr["af_mid"] is not None and fr["an_mid"] else fr["af"]
    haps: dict[str, list[list[int]]] = {}

    def gamete(par):
        h, out = rng.randrange(2), []
        for k in range(len(sites)):
            if k and rng.random() < recomb:
                h = 1 - h
            out.append(par[h][k])
        return out

    def get(s):
        if s not in haps:
            f, m = PARENTS[s]
            haps[s] = ([[int(rng.random() < af(fr)) for *_, fr in sites] for _ in range(2)]
                       if f is None else [gamete(get(f)), gamete(get(m))])
        return haps[s]

    out = {s: [list(get(s)[0]), list(get(s)[1])] for s in SEQ}
    for s in SEQ:                                     # variante sembrada recesiva
        g = SEEDED_GT[s]
        out[s][0][seed_idx] = int(g >= 1)
        out[s][1][seed_idx] = int(g == 2)
    return out                                        # haplotipos por individuo


def write_bams(work: Path, chrom: str, length: int, seq: str, start: int, sites, haps,
               depth: int, rng: random.Random) -> dict[str, Path]:
    """Cada lectura sale de un haplotipo y lleva sus alelos en TODOS los sitios que cubre."""
    import bisect
    import pysam
    positions = [p for p, *_ in sites]
    bams = {}
    for s in SEQ:
        hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": chrom, "LN": length}],
               "RG": [{"ID": "rg", "SM": s}]}
        reads, n = [], 0
        for k, (p, _r, alt, _fr) in enumerate(sites):
            for _ in range(depth):
                rs = max(start, p - rng.randint(0, 99))
                if rs + 100 > start + len(seq):
                    rs = start + len(seq) - 100
                r = list(seq[rs - start: rs - start + 100])
                h = haps[s][rng.randrange(2)]
                lo, hi = bisect.bisect_left(positions, rs), bisect.bisect_right(positions, rs + 99)
                for j in range(lo, hi):
                    if h[j]:
                        r[positions[j] - rs] = sites[j][2]
                a = pysam.AlignedSegment()
                a.query_name = f"{s}_{n}"; n += 1
                a.query_sequence = "".join(r)
                a.flag = 0; a.reference_id = 0; a.reference_start = rs - 1
                a.mapping_quality = 60; a.cigar = [(0, 100)]
                a.query_qualities = pysam.qualitystring_to_array("I" * 100)
                a.set_tag("RG", "rg")
                reads.append(a)
        reads.sort(key=lambda x: x.reference_start)
        path = work / "bams" / f"{s}.bam"
        path.parent.mkdir(parents=True, exist_ok=True)
        with pysam.AlignmentFile(str(path), "wb", header=hdr) as o:
            for a in reads:
                o.write(a)
        pysam.index(str(path))
        bams[s] = path
    return bams


def main(argv: list[str] | None = None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--workdir", type=Path, default=Path.home() / "pga_e2")
    ap.add_argument("--chrom", default="chr22")
    ap.add_argument("--start", type=int, default=42123940)
    ap.add_argument("--end", type=int, default=42133941)
    ap.add_argument("--label", default="CYP2D6")
    ap.add_argument("--min-site-af", type=float, default=0.005,
                    help="SNV de gnomAD a simular (además de la sembrada)")
    ap.add_argument("--depth", type=int, default=30)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--max-af", type=float, default=0.01)
    ap.add_argument("--seq-fasta", type=Path, help="FASTA de la ROI (en vez de UCSC)")
    ap.add_argument("--gnomad-json", type=Path, help="variantes gnomAD (en vez de la API)")
    a = ap.parse_args(argv)

    from publicgenomicagent.knowledge import gnomad as gn
    from publicgenomicagent.knowledge.reference import fetch_sequence

    work = a.workdir.expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(a.seed)
    length = GRCH38_LEN.get(a.chrom)
    if length is None:
        raise SystemExit(f"longitud GRCh38 de {a.chrom} no registrada")

    log("1/6 secuencia GRCh38 de la ROI")
    seq = read_slice(a.seq_fasta or fetch_sequence(a.chrom, a.start, a.end))
    if len(seq) != a.end - a.start + 1:
        raise SystemExit(f"secuencia de {len(seq)} pb para una ROI de {a.end - a.start + 1}")
    ref = build_reference(work, a.chrom, a.start, seq, length)

    log("2/6 gnomAD v4 y comprobación de REF")
    if a.gnomad_json:
        variants = json.loads(a.gnomad_json.read_text())
        gn.AF_CACHE.mkdir(parents=True, exist_ok=True)   # misma caché que usará la tool
        (gn.AF_CACHE / f"{a.chrom.removeprefix('chr')}_{a.start}_{a.end}.json").write_text(
            json.dumps(variants))
    else:
        variants = gn.fetch_region_af(a.chrom, a.start, a.end)
    snvs, mism = snvs_matching_reference(variants, seq, a.start)
    log(f"   {len(snvs)} SNV con REF coincidente, {mism} discordantes")
    if mism > max(5, 0.01 * (len(snvs) + mism)):
        raise SystemExit("demasiados REF discordantes: referencia o coordenadas incorrectas")

    rare = [x for x in snvs if 0 < (x[3]["af"] or 0) < 1e-4
            and (x[3]["af_mid"] or 0) < 1e-3 and 150 < x[0] - a.start < len(seq) - 150]
    if not rare:
        raise SystemExit("no hay SNV raras en gnomAD para sembrar")
    seeded = rare[rng.randrange(len(rare))]
    chosen = {x[0]: x for x in snvs if max(x[3]["af"] or 0, x[3]["af_mid"] or 0) >= a.min_site_af
              and 150 < x[0] - a.start < len(seq) - 150}
    chosen[seeded[0]] = seeded
    sites = [chosen[p] for p in sorted(chosen)]
    seed_idx = sorted(chosen).index(seeded[0])
    log(f"3/6 familia simulada sobre {len(sites)} SNV reales; sembrada "
        f"{a.chrom}:{seeded[0]} {seeded[1]}>{seeded[2]} (AF {seeded[3]['af']:.2e})")
    haps = simulate(sites, seed_idx, rng)

    log("4/6 BAMs")
    bams = write_bams(work, a.chrom, length, seq, a.start, sites, haps, a.depth, rng)
    (work / "family.ped").write_text(PED)
    truth = {"chrom": a.chrom, "pos": seeded[0], "ref": seeded[1], "alt": seeded[2],
             "gnomad_af": seeded[3]["af"], "gnomad_af_mid": seeded[3]["af_mid"],
             "n_sites": len(sites), "n_common_sites": sum(
                 max(x[3]["af"] or 0, x[3]["af_mid"] or 0) > a.max_af for x in sites)}
    (work / "truth.json").write_text(json.dumps(truth, indent=1))

    log("5/6 pga run-family --max-af")
    import os
    from publicgenomicagent.cli.main import app
    os.chdir(work)
    args = ["run-family", "--ped", str(work / "family.ped"), "--reference", str(ref),
            "--region", f"{a.chrom}:{a.start}-{a.end}", "--label", a.label,
            "--case-id", "e2", "--proband", "P1", "--max-af", str(a.max_af)]
    for s, b in bams.items():
        args += ["--bam", f"{s}={b}"]
    try:
        app(args, standalone_mode=False)
    except SystemExit as e:
        if e.code:
            raise

    log("6/6 evaluación")
    st = json.loads((work / "results" / "e2.session.json").read_text())
    calls = {c["tool_name"]: c for c in st["tool_calls"]}
    mf = calls["mendelian_filter"]["output"]
    rep = json.loads(Path(mf["report_json"]).read_text())
    rec_vcf = Path(mf["auto_rec_hom_vcf"])
    import gzip
    rec = [l.split("\t")[:5] for l in gzip.open(rec_vcf, "rt") if not l.startswith("#")]
    in_rec = any(int(p) == seeded[0] and r == seeded[1] and x == seeded[2] for _, p, _, r, x in rec)
    ranking = (calls.get("prioritize_variants") or {}).get("output", {}).get("ranking", [])
    rank = next((v["rank"] for v in ranking if v["pos"] == seeded[0]), None)
    freq = calls["annotate_frequency"]["output"]
    result = {
        "seeded_in_recessive_candidates": in_rec,
        "n_recessive_candidates": rep["counts"]["auto_rec_hom"],
        "seeded_rank_prioritized": rank,
        "rarity_filter": rep["rarity_filter"],
        "annotate_frequency": freq["message"],
        "gnomad_method": freq["method"],
        "n_ref_mismatch": freq.get("n_ref_mismatch"),
        "comp_het": rep["counts"]["comp_het"],
    }
    passed = (in_rec and rep["rarity_filter"]["removed_common"] > 0
              and rep["counts"]["auto_rec_hom"] <= 5 and not freq.get("n_ref_mismatch"))
    result["E2_passed"] = passed
    md = ["# E2: caso realista (GRCh38 real, variante recesiva sembrada)", "",
          f"- ROI: {a.chrom}:{a.start}-{a.end} ({a.label}); referencia GRCh38 real (UCSC) en la ROI",
          f"- SNV de gnomAD simuladas: {len(sites)} (de ellas comunes > {a.max_af}: "
          f"{truth['n_common_sites']}); REF discordantes con gnomAD: {mism}",
          f"- Variante sembrada: {a.chrom}:{seeded[0]} {seeded[1]}>{seeded[2]}, "
          f"AF gnomAD {seeded[3]['af']:.2e}",
          f"- Filtro de rareza: {rep['rarity_filter']}",
          f"- Candidatas recesivas homocigotas: {rep['counts']['auto_rec_hom']} "
          f"(sembrada incluida: {'sí' if in_rec else 'NO'}); posición en la priorización: {rank}",
          f"- Heterocigoto compuesto: {rep['counts']['comp_het']} en trans, "
          f"{rep['counts']['comp_het_unphased']} sin fase",
          f"- **E2 {'superado' if passed else 'NO superado'}**", ""]
    (work / "e2_report.md").write_text("\n".join(md))
    print("\n".join(md))
    (work / "e2_result.json").write_text(json.dumps(result, indent=1, ensure_ascii=False))
    return result


if __name__ == "__main__":
    sys.exit(0 if main().get("E2_passed") else 1)
