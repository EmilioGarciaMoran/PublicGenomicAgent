#!/usr/bin/env python3
"""Benchmark justo: alineamiento lineal frente a pangenoma local.

Las dos ramas usan las MISMAS lecturas, el MISMO alineador (vg giraffe)
y el MISMO llamador (bcftools mpileup + call sobre el BAM proyectado a
la referencia lineal). Solo cambia el grafo contra el que se alinea:

- lineal:    grafo construido solo con la referencia;
- pangenoma: referencia + alelos de la cohorte (VCF). Las variantes
             objetivo NUNCA están en el VCF de la cohorte: el grafo no
             puede "saber" la respuesta.

Escenarios sobre una referencia sintética (semilla fija), muestra
homocigota:

1. control:     SNV nueva en región única (ambas ramas deben verla);
2. divergente:  SNV nueva dentro de un haplotipo muy distinto de la
                referencia (SNV cada ~12 pb y pequeñas indels, todas en
                la cohorte): sesgo de referencia;
3. duplicación: SNV nueva en la copia 1 de una duplicación segmentaria
                (copias con ~1 % de diferencias, PSV); la muestra lleva
                en la copia 1 un tramo convertido con las bases de la
                copia 2 (en la cohorte): las lecturas parecen de la copia 2.

Salida: tabla Markdown (y JSON con --json) con, por rama y escenario, si
se llama la variante objetivo, genotipo, profundidad y MAPQ medio en la
posición, y falsos positivos en todo el contig.

Uso:  python scripts/benchmark_local_pangenome.py --workdir DIR [--json out.json]
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import statistics
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from publicgenomicagent.env.micromamba import bin_path  # noqa: E402

CHROM = "chrS"
L = 30_000
SEED = 7
BASES = "ACGT"


def other(b: str, rng: random.Random) -> str:
    return rng.choice([x for x in BASES if x != b])


def design(rng: random.Random, spacing: int = 12):
    """Referencia, alelos de cohorte (que la muestra lleva) y objetivos."""
    ref = [rng.choice(BASES) for _ in range(L)]
    # duplicación: copia 1 en [4000, 7000), copia 2 en [20000, 23000), ~1 % PSV
    c1, c2, dl = 4000, 20000, 3000
    for i in range(dl):
        ref[c2 + i] = ref[c1 + i]
    psv = sorted(rng.sample(range(dl), 30))
    for i in psv:
        ref[c2 + i] = other(ref[c1 + i], rng)

    cohort: dict[int, tuple[str, str]] = {}           # pos0 -> (ref, alt), alelos de la muestra
    # haplotipo divergente en [12000, 13500): SNV cada ~12 pb + 4 indels pequeñas
    for p in range(12000, 13500, spacing):
        q = p + rng.randrange(0, max(1, spacing // 2))
        cohort[q] = (ref[q], other(ref[q], rng))
    for p in (12210, 12630, 13050, 13390):
        for d in range(-3, 9):                      # sin SNV pegadas a la indel
            cohort.pop(p + d, None)
        if p in (12210, 13050):                     # deleción de 6 pb
            cohort[p] = ("".join(ref[p:p + 7]), ref[p])
        else:                                       # inserción de 5 pb
            cohort[p] = (ref[p], ref[p] + "".join(rng.choice(BASES) for _ in range(5)))
    # conversión génica: en la copia 1, PSV dentro de ±400 pb del objetivo -> base de la copia 2
    t_dup = c1 + 1500
    for i in psv:
        if abs(c1 + i - t_dup) <= 400:
            cohort[c1 + i] = (ref[c1 + i], ref[c2 + i])

    targets = {"control": 27000, "divergent": 12750, "duplication": t_dup}
    for t in targets.values():
        for d in range(-2, 3):
            cohort.pop(t + d, None)
    alts = {k: (t, ref[t], other(ref[t], rng)) for k, t in targets.items()}
    return "".join(ref), cohort, alts, {"copy1": (c1, c1 + dl), "copy2": (c2, c2 + dl), "psv": len(psv)}


def sample_seq(ref: str, cohort, alts) -> str:
    edits = {p: (r, a) for p, (r, a) in cohort.items()}
    for t, r, a in alts.values():
        edits[t] = (r, a)
    out, i = [], 0
    for p in sorted(edits):
        r, a = edits[p]
        out.append(ref[i:p]); out.append(a); i = p + len(r)
    out.append(ref[i:])
    return "".join(out)


def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], check=True, capture_output=True, text=True, **kw)


def write_fa(path: Path, name: str, seq: str):
    with path.open("w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")


def write_vcf(path: Path, cohort):
    lines = ["##fileformat=VCFv4.2", f"##contig=<ID={CHROM},length={L}>",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for p in sorted(cohort):
        r, a = cohort[p]
        lines.append(f"{CHROM}\t{p + 1}\t.\t{r}\t{a}\t.\tPASS\t.")
    path.write_text("\n".join(lines) + "\n")


def calls(vcf: Path):
    out = {}
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            fmt = dict(zip(c[8].split(":"), c[9].split(":")))
            out[(int(c[1]), c[3], c[4])] = (fmt.get("GT"), float(c[5]))
    return out


def mapq_depth(bam: Path, pos1: int, samtools: Path):
    out = run([samtools, "view", bam, f"{CHROM}:{pos1}-{pos1}"]).stdout.splitlines()
    mq = [int(l.split("\t")[4]) for l in out]
    return len(mq), (statistics.mean(mq) if mq else 0.0), sum(m == 0 for m in mq)


def branch(name: str, wd: Path, ref_fa: Path, vcf_gz: Path | None, r1: Path, r2: Path):
    vg = bin_path("pga-pangenome", "vg")
    samtools, bcftools = bin_path("pga-hts", "samtools"), bin_path("pga-hts", "bcftools")
    pre = wd / name
    idx = [vg, "autoindex", "--workflow", "giraffe", "-r", ref_fa, "-p", pre, "-T", wd / "tmp"]
    if vcf_gz:
        idx += ["-v", vcf_gz]
    (wd / "tmp").mkdir(exist_ok=True)
    run(idx)
    raw = wd / f"{name}.raw.bam"
    with raw.open("wb") as f:
        subprocess.run([str(vg), "giraffe", "-Z", f"{pre}.giraffe.gbz", "-m", f"{pre}.shortread.withzip.min",
                        "-z", f"{pre}.shortread.zipcodes", "-d", f"{pre}.dist", "-f", str(r1), "-f", str(r2), "-o", "BAM",
                        "--sample", "S", "--read-group", "S"],
                       check=True, stdout=f, stderr=subprocess.PIPE)
    bam = wd / f"{name}.bam"
    run([samtools, "sort", "-o", bam, raw]); run([samtools, "index", bam])
    vcf = wd / f"{name}.vcf.gz"
    out = {}
    for q, v in ((0, vcf), (20, wd / f"{name}.mq20.vcf.gz")):
        mp = subprocess.Popen([str(bcftools), "mpileup", "-f", str(ref_fa), "-a", "AD,DP", "-q", str(q),
                               "-d", "10000", str(bam)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        run([bcftools, "call", "-mv", "-Oz", "-o", v], stdin=mp.stdout)
        mp.wait()
        out[q] = v
    return bam, out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--coverage", type=int, default=30)
    ap.add_argument("--divergent-spacing", type=int, default=12,
                    help="pb entre SNV de cohorte en el haplotipo divergente (12 ~ 8 %%, 100 ~ 1 %%)")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--json")
    ap.add_argument("--sweep", action="store_true",
                    help="10 semillas x divergencias 100/50/25/12 pb; imprime la tabla agregada")
    a = ap.parse_args(argv)
    if a.sweep:
        return sweep(Path(a.workdir), a.coverage, a.json)
    wd = Path(a.workdir).expanduser().resolve(); wd.mkdir(parents=True, exist_ok=True)
    rng = random.Random(a.seed)
    ref, cohort, alts, info = design(rng, a.divergent_spacing)
    samp = sample_seq(ref, cohort, alts)

    ref_fa, samp_fa = wd / "ref.fa", wd / "sample.fa"
    write_fa(ref_fa, CHROM, ref); write_fa(samp_fa, CHROM, samp)
    samtools, bcftools = bin_path("pga-hts", "samtools"), bin_path("pga-hts", "bcftools")
    run([samtools, "faidx", ref_fa])
    vcf = wd / "cohort.vcf"; write_vcf(vcf, cohort)
    vcf_gz = wd / "cohort.vcf.gz"
    run([bcftools, "view", vcf, "-Oz", "-o", vcf_gz]); run([bcftools, "index", "-t", vcf_gz])
    # lecturas pareadas 2x150, inserto 400, error 0,2 %, sin mutaciones añadidas por wgsim
    n = a.coverage * len(samp) // 300
    r1, r2 = wd / "r1.fq", wd / "r2.fq"
    run([bin_path("pga-hts", "wgsim"), "-S", a.seed, "-N", n, "-1", 150, "-2", 150, "-d", 400,
         "-s", 40, "-e", 0.002, "-r", 0, "-R", 0, samp_fa, r1, r2])

    truth = {(p + 1, r, al) for p, (r, al) in cohort.items()} | {(t + 1, r, al) for t, r, al in alts.values()}
    res = {"design": {**info, "cohort_variants": len(cohort), "coverage": a.coverage,
                      "reads_pairs": n, "seed": a.seed, "targets_in_cohort_vcf": 0,
                      "divergent_spacing_bp": a.divergent_spacing}, "branches": {}}
    for name, v in (("linear", None), ("local_pangenome", vcf_gz)):
        bam, out = branch(name, wd, ref_fa, v, r1, r2)
        cs, cs20 = calls(out[0]), calls(out[20])
        row = {}
        for k, (t, r, al) in alts.items():
            d, mq, mq0 = mapq_depth(bam, t + 1, samtools)
            gt, q = cs.get((t + 1, r, al), (None, None))
            row[k] = {"pos": t + 1, "called": gt is not None, "gt": gt, "qual": q,
                      "called_mapq20": (t + 1, r, al) in cs20,
                      "depth": d, "mean_mapq": round(mq, 1), "mapq0": mq0}
        row["false_positives"] = sum(1 for key in cs if key not in truth)
        row["cohort_alleles_recalled"] = sum(1 for key in cs if key in truth) - sum(r["called"] for r in row.values() if isinstance(r, dict))
        res["branches"][name] = row

    print(f"# Linear vs local pangenome ({CHROM}, {L} bp, seed {a.seed}, {a.coverage}x)\n")
    print(f"Cohort alleles in graph: {len(cohort)}; target variants in cohort VCF: 0\n")
    print(f"Divergent haplotype: one cohort SNV every ~{a.divergent_spacing} bp\n")
    print("| Scenario | Branch | Called | Called (MAPQ>=20) | GT | Depth | Mean MAPQ |")
    print("|---|---|---|---|---|---|---|")
    for k in alts:
        for b in res["branches"]:
            r = res["branches"][b][k]
            print(f"| {k} | {b} | {'yes' if r['called'] else 'no'} | "
                  f"{'yes' if r['called_mapq20'] else 'no'} | {r['gt'] or '—'} | "
                  f"{r['depth']} | {r['mean_mapq']} |")
    for b, r in res["branches"].items():
        print(f"\n{b}: false positives {r['false_positives']}, cohort alleles recalled "
              f"{r['cohort_alleles_recalled']}/{len(cohort)}")
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=2))
    return 0


def sweep(workdir: Path, coverage: int, json_out: str | None,
          spacings=(100, 50, 25, 12), seeds=range(1, 11)) -> int:
    import contextlib
    import io
    import shutil
    agg: dict = {}
    for sp in spacings:
        for sd in seeds:
            wd = workdir / f"sp{sp}_seed{sd}"
            j = workdir / f"sp{sp}_seed{sd}.json"
            with contextlib.redirect_stdout(io.StringIO()):
                main(["--workdir", str(wd), "--coverage", str(coverage), "--divergent-spacing", str(sp),
                      "--seed", str(sd), "--json", str(j)])
            shutil.rmtree(wd, ignore_errors=True)
            r = json.loads(j.read_text())
            for b, row in r["branches"].items():
                for k in ("control", "divergent", "duplication"):
                    agg.setdefault((k, sp, b), []).append(row[k])
                agg.setdefault(("fp", sp, b), []).append(row["false_positives"])
    n = len(list(seeds))
    print(f"# Linear vs local pangenome: {n} seeds per divergence, {coverage}x\n")
    print("| Scenario | Divergent SNV every | Branch | Called | Called (MAPQ>=20) | Mean depth | Mean MAPQ |")
    print("|---|---|---|---|---|---|---|")
    table = []
    for k in ("control", "divergent", "duplication"):
        for sp in spacings:
            for b in ("linear", "local_pangenome"):
                v = agg[(k, sp, b)]
                row = {"scenario": k, "spacing": sp, "branch": b, "n": len(v),
                       "called": sum(x["called"] for x in v),
                       "called_mapq20": sum(x["called_mapq20"] for x in v),
                       "mean_depth": round(statistics.mean(x["depth"] for x in v), 1),
                       "mean_mapq": round(statistics.mean(x["mean_mapq"] for x in v), 1)}
                table.append(row)
                print(f"| {k} | {sp} bp | {b} | {row['called']}/{n} | {row['called_mapq20']}/{n} | "
                      f"{row['mean_depth']} | {row['mean_mapq']} |")
    fps = {sp: {b: sum(agg[("fp", sp, b)]) for b in ("linear", "local_pangenome")} for sp in spacings}
    print("\nFalse positives (sum over seeds): " + "; ".join(
        f"{sp} bp: linear {v['linear']}, pangenome {v['local_pangenome']}" for sp, v in fps.items()))
    if json_out:
        Path(json_out).write_text(json.dumps({"table": table, "false_positives": fps}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
