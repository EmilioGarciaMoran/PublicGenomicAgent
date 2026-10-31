#!/usr/bin/env python3
"""Calibra la verificación familiar (IBS local) con el trío GIAB recortado.

Datos (sin descargar BAMs completos: lectura remota por región vía índice):
  - Exomas Illumina de Oslo University Hospital (Agilent SureSelect v5,
    bwa-mem, GRCh37) de HG002 (hijo), HG003 (padre), HG004 (madre) del
    trío askenazí y HG005 (trío chino) como "padre falso" de otra familia.
    URLs y MD5 tomados de github.com/genome-in-a-bottle/giab_data_indexes.
  - Referencia hg19 por cromosoma (UCSC), renombrada a los contigs del BAM.
  - SNV comunes de 1000 Genomes fase 3 (los GIAB no están en 1000G).

Pasos (cada uno se salta si su salida ya existe en --workdir):
  1. cabeceras remotas: nombres de contig, longitudes y SM;
  2. referencia de los cromosomas pedidos, con comprobación de longitudes;
  3. sub-BAM por muestra de esos cromosomas, SM reetiquetado al ID GIAB;
  4. regiones cubiertas en HG002 (profundidad >= --min-depth);
  5. panel: SNV bialélicos de 1000G con AF en [--af-min, --af-max] dentro
     de las regiones cubiertas (VCF de 1000G descargado una vez, con
     reintentos y reanudación, y filtrado en local);
  6. `verify_relatedness` real sobre las cuatro muestras;
  7. submuestreo de paneles de --ks sitios (--reps repeticiones) en dos
     escenarios: trío verdadero (debe pasar) y HG005 declarado padre de
     HG002 (debe fallar);
  8. informe en --workdir/report.md y report.json.

Uso típico (en la máquina con red y pga-hts):
  pga-env
  python scripts/calibrate_relatedness_giab.py --workdir ~/pga_giab_calib

Limitación conocida: el "padre falso" es de otra ascendencia (HG005, chino);
detectar un no-progenitor de la misma población es algo más difícil. No hay
en GIAB una mujer no emparentada con exoma de Oslo para probar "madre falsa".
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

from publicgenomicagent.env.micromamba import bin_path
from publicgenomicagent.tools.base import VerifyRelatednessInput
from publicgenomicagent.tools.pedigree import parse_ped_text
from publicgenomicagent.tools.relatedness import load_sites, parse_genotypes, verify_relatedness
from publicgenomicagent.tools.relatedness_calibration import (
    Scenario,
    calibrate,
    markdown_report,
    pairwise_summary,
)

GIAB_OSLO_EXOMES = {
    "HG002": {
        "bam": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam",
        "bam_md5": "c80f0cab24bfaa504393457b8f7191fa",
        "bai": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai",
        "bai_md5": "d4fea426c3e2e9a71bb92e6526b4df6f",
    },
    "HG003": {
        "bam": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam",
        "bam_md5": "79ce0edc9363710030e973bd12af5423",
        "bai": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bai",
        "bai_md5": "0512b54f2b66f3a04653be9f5975ae7",  # 31 caracteres en el índice oficial
    },
    "HG004": {
        "bam": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam",
        "bam_md5": "39e90d925fff5ea7dd6dab255f299581",
        "bai": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bai",
        "bai_md5": "8914bfb6fa6bd304192f2c9e13e903f4",
    },
    "HG005": {
        "bam": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG005-EEogPU_v02-KIT-Av5_CGCATACA_L008.posiSrt.markDup.bam",
        "bam_md5": "eceb2ace6aa97dcd48628ad7a3e5698e",
        "bai": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG005-EEogPU_v02-KIT-Av5_CGCATACA_L008.posiSrt.markDup.bai",
        "bai_md5": "a864f73b7598145df835297c019da5f1",
    },
}

UCSC_HG19 = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr{chrom}.fa.gz"
ONEKG_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"

PED_ALL = "GIAB HG003 0 0 1 1\nGIAB HG004 0 0 2 1\nGIAB HG002 HG003 HG004 1 2\nGIAB HG005 0 0 1 1\n"
PED_TRUE = "GIAB HG003 0 0 1 1\nGIAB HG004 0 0 2 1\nGIAB HG002 HG003 HG004 1 2\n"
PED_FALSE_FATHER = "GIAB HG005 0 0 1 1\nGIAB HG004 0 0 2 1\nGIAB HG002 HG005 HG004 1 2\n"


def log(msg: str) -> None:
    print(f"[calib] {msg}", flush=True)


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


def fetch(src: str, dest: Path, retries: int = 5, check=None) -> Path:
    """Descarga (http/https/ftp) o copia (ruta local) `src` a `dest`.

    Robusto frente a cortes: reintenta y, en HTTP, retoma desde lo ya
    descargado (cabecera Range) en `dest.part`. `check(path)` valida el
    fichero completo (p.ej. integridad BGZF); si falla se borra y reintenta.
    """
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    if not re.match(r"^(https?|ftp)://", src):
        shutil.copyfile(src, tmp)
        tmp.rename(dest)
        return dest
    last: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            have = tmp.stat().st_size if tmp.exists() else 0
            req = urllib.request.Request(src, headers={"Range": f"bytes={have}-"} if have else {})
            with urllib.request.urlopen(req, timeout=300) as r:
                resumed = have and getattr(r, "status", 200) == 206
                with tmp.open("ab" if resumed else "wb") as f:
                    shutil.copyfileobj(r, f, length=1 << 20)
                total = r.headers.get("Content-Range", "").rsplit("/", 1)[-1] if resumed \
                    else r.headers.get("Content-Length")
            if total and total.isdigit() and tmp.stat().st_size != int(total):
                raise OSError(f"descarga incompleta: {tmp.stat().st_size}/{total} bytes")
            if check is not None and not check(tmp):
                tmp.unlink()
                raise OSError("el fichero descargado no supera la comprobación de integridad")
            tmp.rename(dest)
            return dest
        except Exception as e:  # noqa: BLE001
            last = e
            log(f"   descarga de {Path(src).name}: intento {attempt}/{retries} falló ({e}); reintento")
    raise SystemExit(f"no se pudo descargar {src}: {last}")


def bam_header(samtools: str, src: dict) -> dict:
    out = run([samtools, "view", "-H", "-X", src["bam"], src["bai"]],
              capture_output=True, text=True).stdout
    sq, sms = {}, []
    for line in out.splitlines():
        f = line.split("\t")
        if f[0] == "@SQ":
            d = dict(x.split(":", 1) for x in f[1:] if ":" in x)
            sq[d["SN"]] = int(d["LN"])
        elif f[0] == "@RG":
            sms += [x[3:] for x in f[1:] if x.startswith("SM:")]
    return {"sq": sq, "sm": sorted(set(sms)), "text": out}


def contig(chrom: str, chr_prefix: bool) -> str:
    c = chrom.removeprefix("chr")
    return f"chr{c}" if chr_prefix else c


def build_reference(work: Path, chroms: list[str], prefix: bool, sq: dict, tmpl: str,
                    samtools: str) -> Path:
    ref = work / "ref" / "grch37_subset.fa"
    if ref.exists() and Path(str(ref) + ".fai").exists():
        return ref
    ref.parent.mkdir(parents=True, exist_ok=True)
    with ref.open("w") as out:
        for c in chroms:
            gz = fetch(tmpl.format(chrom=c), work / "ref" / f"chr{c}.fa.gz")
            name = contig(c, prefix)
            with gzip.open(gz, "rt") as f:
                for line in f:
                    out.write(f">{name}\n" if line.startswith(">") else line.upper())
    run([samtools, "faidx", ref])
    lengths = {l.split("\t")[0]: int(l.split("\t")[1]) for l in Path(str(ref) + ".fai").read_text().splitlines()}
    for c in chroms:
        name = contig(c, prefix)
        if name not in sq or sq[name] != lengths.get(name):
            ref.unlink()
            raise SystemExit(f"longitud de {name}: BAM {sq.get(name)} vs referencia {lengths.get(name)} "
                             "(¿otra build?). Se aborta.")
    return ref


def extract_subbam(work: Path, sid: str, src: dict, regions: list[str], hdr: dict,
                   samtools: str) -> Path:
    final = work / "bams" / f"{sid}.bam"
    if final.exists() and Path(str(final) + ".bai").exists():
        return final
    final.parent.mkdir(parents=True, exist_ok=True)
    raw = work / "bams" / f"{sid}.raw.bam"
    log(f"{sid}: lectura remota de {', '.join(regions)} …")
    run([samtools, "view", "-b", "-X", "-o", raw, src["bam"], src["bai"], *regions])
    new_hdr = "\n".join(
        re.sub(r"\tSM:[^\t]+", f"\tSM:{sid}", l) if l.startswith("@RG") else l
        for l in hdr["text"].splitlines()) + "\n"
    hfile = work / "bams" / f"{sid}.header.sam"
    hfile.write_text(new_hdr)
    with final.open("wb") as f:
        run([samtools, "reheader", hfile, raw], stdout=f)
    raw.unlink()
    run([samtools, "index", final])
    log(f"{sid}: sub-BAM {final.stat().st_size / 1e6:.1f} MB")
    return final


def covered_regions(work: Path, bam: Path, min_depth: int, samtools: str) -> Path:
    """Intervalos con profundidad >= min_depth (1-based, CHROM BEG END)."""
    out = work / "panel" / "covered.tsv"
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    p = subprocess.Popen([samtools, "depth", "-q", "20", "-Q", "20", str(bam)],
                         stdout=subprocess.PIPE, text=True)
    n = 0
    with out.open("w") as f:
        cur = None
        for line in p.stdout:
            c, pos, d = line.split("\t")
            pos, d = int(pos), int(d)
            if d < min_depth:
                continue
            if cur and cur[0] == c and pos == cur[2] + 1:
                cur[2] = pos
            else:
                if cur:
                    f.write(f"{cur[0]}\t{cur[1]}\t{cur[2]}\n"); n += 1
                cur = [c, pos, pos]
        if cur:
            f.write(f"{cur[0]}\t{cur[1]}\t{cur[2]}\n"); n += 1
    if p.wait() != 0:
        raise SystemExit("samtools depth falló")
    log(f"regiones cubiertas (>= {min_depth}x en HG002): {n}")
    return out


def onekg_vcf(base: str, chrom: str) -> str:
    """Localiza el VCF por cromosoma de 1000G fase 3 (el sufijo v5a/v5b varía)."""
    c = chrom.removeprefix("chr")
    pat = re.compile(rf"ALL\.chr{re.escape(c)}\.phase3_[A-Za-z0-9_.]+?\.genotypes\.vcf\.gz(?!\.tbi)")
    if re.match(r"^(https?|ftp)://", base):
        with urllib.request.urlopen(base, timeout=120) as r:
            listing = r.read().decode("utf-8", "replace")
        names = sorted(set(m.group(0) for m in pat.finditer(listing)))
    else:
        names = sorted(p.name for p in Path(base).iterdir() if pat.fullmatch(p.name))
    if not names:
        raise SystemExit(f"no encuentro el VCF de 1000G para chr{c} en {base}")
    return base.rstrip("/") + "/" + names[-1]


def build_sites(work: Path, chroms: list[str], prefix: bool, covered: Path, base: str,
                af_min: float, af_max: float, bcftools: str) -> Path:
    sites = work / "panel" / "sites.tsv"
    if sites.exists():
        return sites
    regions_1kg = work / "panel" / "covered_1kg.tsv"
    regions_1kg.write_text("".join(
        f"{l.split(chr(9))[0].removeprefix('chr')}\t" + "\t".join(l.split("\t")[1:])
        for l in covered.read_text().splitlines(keepends=True)))
    seen, rows = set(), []
    for c in chroms:
        url = onekg_vcf(base, c)
        name = url.rsplit("/", 1)[-1]
        log(f"1000G chr{c}: {name} (descarga local una vez; ~200 MB)")
        # Se descarga entero: consultar en remoto miles de regiones pequeñas
        # obliga a miles de saltos HTTP y un solo corte rompe la lectura.
        local = fetch(url, work / "onekg" / name, check=lambda p: bgzip_ok(p, bcftools))
        fetch(url + ".tbi", work / "onekg" / (name + ".tbi"))
        view = subprocess.Popen(
            [bcftools, "view", "-R", str(regions_1kg), "-G", "-v", "snps", "-m2", "-M2",
             "-i", f"INFO/AF>={af_min} && INFO/AF<={af_max}", "-Ou", str(local)],
            stdout=subprocess.PIPE)
        q = subprocess.run([bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n"],
                           stdin=view.stdout, capture_output=True, text=True, check=True)
        if view.wait() != 0:
            raise SystemExit(f"bcftools view falló en {url}")
        for line in q.stdout.splitlines():
            ch, pos, ref, alt = line.split("\t")
            key = (contig(ch, prefix), int(pos))
            if key not in seen:
                seen.add(key)
                rows.append(f"{key[0]}\t{key[1]}\t{ref}\t{alt}\n")
    sites.write_text("".join(rows))
    log(f"panel de sitios: {len(rows)} SNV comunes en regiones cubiertas")
    return sites


def bgzip_ok(path: Path, bcftools: str) -> bool:
    """Integridad BGZF completa (bgzip -t, del mismo entorno que bcftools)."""
    bgzip = str(Path(bcftools).with_name("bgzip"))
    return subprocess.run([bgzip, "-t", str(path)], capture_output=True).returncode == 0


def tool_version(binary: str) -> str:
    out = subprocess.run([binary, "--version"], capture_output=True, text=True).stdout
    return " / ".join(out.splitlines()[:2])


def main(argv: list[str] | None = None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--workdir", type=Path, default=Path.home() / "pga_giab_calibration")
    ap.add_argument("--chroms", default="21,22")
    ap.add_argument("--min-depth", type=int, default=20)
    ap.add_argument("--af-min", type=float, default=0.2)
    ap.add_argument("--af-max", type=float, default=0.8)
    ap.add_argument("--ks", default="50,100,150,200,300,500,1000")
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sources", type=Path, help="JSON {ID: {bam, bai}} en lugar de GIAB")
    ap.add_argument("--ref-template", default=UCSC_HG19)
    ap.add_argument("--onekg-base", default=ONEKG_BASE)
    a = ap.parse_args(argv)

    samtools = str(bin_path("pga-hts", "samtools"))
    bcftools = str(bin_path("pga-hts", "bcftools"))
    for b in (samtools, bcftools):
        if not Path(b).exists():
            raise SystemExit(f"no encuentro {b}: activa pga-env / crea pga-hts")
    sources = json.loads(a.sources.read_text()) if a.sources else GIAB_OSLO_EXOMES
    chroms = [c.strip().removeprefix("chr") for c in a.chroms.split(",") if c.strip()]
    work = a.workdir.expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)

    log("1/8 cabeceras remotas")
    headers = {sid: bam_header(samtools, src) for sid, src in sources.items()}
    first = headers["HG002"]
    prefix = any(n.startswith("chr") for n in first["sq"])
    for sid, h in headers.items():
        log(f"   {sid}: SM={h['sm']} contigs={'chrN' if prefix else 'N'}")

    log("2/8 referencia")
    ref = build_reference(work, chroms, prefix, first["sq"], a.ref_template, samtools)
    log("3/8 sub-BAMs")
    regions = [contig(c, prefix) for c in chroms]
    bams = {sid: extract_subbam(work, sid, src, regions, headers[sid], samtools)
            for sid, src in sources.items()}
    log("4/8 regiones cubiertas")
    covered = covered_regions(work, bams["HG002"], a.min_depth, samtools)
    log("5/8 panel de sitios (1000G)")
    sites_tsv = build_sites(work, chroms, prefix, covered, a.onekg_base, a.af_min, a.af_max, bcftools)

    log("6/8 verify_relatedness (tool real)")
    _, ped_all = parse_ped_text(PED_ALL)
    vr = verify_relatedness(None, VerifyRelatednessInput(
        bams=[(sid, bams[sid]) for sid in ("HG003", "HG004", "HG002", "HG005")],
        reference_fasta=ref, sites=sites_tsv, pedigree=ped_all, output_dir=work / "verify"))
    log(f"   estado con todos los sitios: {vr.status}")

    log("7/8 submuestreo de paneles")
    sites = load_sites(sites_tsv)
    with gzip.open(vr.genotypes_vcf, "rt") as f:
        samples, gts = parse_genotypes(f, sites, 7)
    ks = [int(k) for k in a.ks.split(",")]
    scenarios = [
        Scenario("trío verdadero (HG003, HG004 -> HG002)", parse_ped_text(PED_TRUE)[1], "passed"),
        Scenario("padre falso (HG005 declarado padre de HG002)", parse_ped_text(PED_FALSE_FATHER)[1], "failed"),
    ]
    cal = calibrate(gts, samples, scenarios, ks, reps=a.reps, seed=a.seed)
    pairs = pairwise_summary(gts, samples)

    log("8/8 informe")
    meta = {
        "fecha": dt.datetime.now().isoformat(timespec="seconds"),
        "datos": "exomas Illumina Oslo (Agilent SureSelect v5, bwa-mem, GRCh37)",
        "cromosomas": ",".join(regions),
        "panel": f"{len(sites)} SNV de 1000G fase 3, AF en [{a.af_min}, {a.af_max}], "
                 f"en regiones >= {a.min_depth}x de HG002",
        "umbrales": "por defecto de verify_relatedness (min_dp 7, min_sites 100, "
                    "max_parent_ibs0 0.01, dup 0.85, deficit_tol 0.2, undeclared 0.3)",
        "repeticiones": f"{a.reps} por tamaño de panel, semilla {a.seed}",
        "samtools": tool_version(samtools),
        "bcftools": tool_version(bcftools),
        "fuentes": "; ".join(f"{s}: {Path(v['bam']).name} (md5 publicado {v.get('bam_md5', '—')})"
                             for s, v in sources.items()),
    }
    report = {"meta": meta, "verify_relatedness_all_sites": vr.model_dump(mode="json"),
              "pairs_all_sites": pairs, "calibration": cal}
    (work / "report.json").write_text(json.dumps(report, indent=1, ensure_ascii=False))
    md = markdown_report(cal, pairs, meta)
    (work / "report.md").write_text(md)
    print(md)
    log(f"informe: {work / 'report.md'}")
    return report


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
