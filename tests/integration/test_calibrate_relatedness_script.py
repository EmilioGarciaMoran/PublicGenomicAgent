"""El script de calibración GIAB, de punta a punta, con fuentes locales sintéticas.

Imita lo que encontrará con los datos reales: contigs "21"/"22" en el BAM y
">chr21" en el FASTA de UCSC, SM con formato de Oslo ("HG002-EEogPU"),
VCF de 1000G con sufijos v5a/v5b, sitios sin cobertura (fuera de "exoma")
y sitios de AF baja que el panel debe descartar.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import random
from pathlib import Path

import pytest

from publicgenomicagent.env.micromamba import bin_path

pysam = pytest.importorskip("pysam")
pytestmark = pytest.mark.integration

SCRIPT = Path(__file__).parents[2] / "scripts" / "calibrate_relatedness_giab.py"
SWAP = {"A": "G", "C": "T", "G": "A", "T": "C"}


def _load_script():
    spec = importlib.util.spec_from_file_location("calib", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _make_sources(root: Path) -> Path:
    rng = random.Random(5)
    chroms = {"21": 90_000, "22": 90_000}
    seqs = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in chroms.items()}
    (root / "ucsc").mkdir(parents=True)
    for c, s in seqs.items():
        with gzip.open(root / "ucsc" / f"chr{c}.fa.gz", "wt") as f:
            f.write(f">chr{c}\n" + "\n".join(s[i:i + 50].lower() for i in range(0, len(s), 50)) + "\n")

    # sitios candidatos (1000G): comunes y raros; ~75% dentro del "exoma"
    sites = []
    for c, s in seqs.items():
        for pos in range(300, chroms[c] - 300, 250):
            af = rng.uniform(0.25, 0.75) if rng.random() < 0.85 else rng.uniform(0.01, 0.1)
            sites.append((c, pos, s[pos - 1], SWAP[s[pos - 1]], af, rng.random() < 0.75))
    onekg = root / "onekg"
    onekg.mkdir()
    for c, suffix in (("21", "v5a"), ("22", "v5b")):
        vcf = onekg / f"ALL.chr{c}.phase3_shapeit2_mvncall_integrated_{suffix}.20130502.genotypes.vcf"
        lines = ["##fileformat=VCFv4.1", f"##contig=<ID={c},length={chroms[c]}>",
                 '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001"]
        for sc, pos, ref, alt, af, _ in sites:
            if sc == c:
                lines.append(f"{c}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\tAF={af:.3f}\tGT\t0|1")
        vcf.write_text("\n".join(lines) + "\n")
        pysam.tabix_index(str(vcf), preset="vcf", force=True)

    # familia: HG003 x HG004 -> HG002; HG005 sin relación
    haps = {}
    for s in ("HG003", "HG004", "HG005"):
        haps[s] = [[int(rng.random() < site[4]) for site in sites] for _ in range(2)]
    haps["HG002"] = [haps["HG003"][rng.randrange(2)], haps["HG004"][rng.randrange(2)]]
    srcs = {}
    for sid, h in haps.items():
        hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
               "SQ": [{"SN": c, "LN": n} for c, n in chroms.items()],
               "RG": [{"ID": "L008", "SM": f"{sid}-EEogPU"}]}
        reads = []
        for k, (c, pos, ref, alt, af, covered) in enumerate(sites):
            if not covered:
                continue
            for j in range(30):
                start = pos - 1 - rng.randint(0, 99)
                r = list(seqs[c][start:start + 100])
                if h[rng.randrange(2)][k]:
                    r[pos - 1 - start] = alt
                a = pysam.AlignedSegment()
                a.query_name = f"{sid}_{k}_{j}"
                a.query_sequence = "".join(r)
                a.flag = 0
                a.reference_id = 0 if c == "21" else 1
                a.reference_start = start
                a.mapping_quality = 60
                a.cigar = [(0, 100)]
                a.query_qualities = pysam.qualitystring_to_array("I" * 100)
                a.set_tag("RG", "L008")
                reads.append(a)
        reads.sort(key=lambda x: (x.reference_id, x.reference_start))
        bam = root / "giab" / f"{sid}.posiSrt.markDup.bam"
        bam.parent.mkdir(exist_ok=True)
        with pysam.AlignmentFile(str(bam), "wb", header=hdr) as o:
            for a in reads:
                o.write(a)
        pysam.index(str(bam))
        bai = bam.with_suffix(".bai")                       # nombre al estilo Oslo
        Path(str(bam) + ".bai").rename(bai)
        srcs[sid] = {"bam": str(bam), "bai": str(bai), "bam_md5": "-"}
    sources = root / "sources.json"
    sources.write_text(json.dumps(srcs))
    return sources


def test_calibration_script_end_to_end(tmp_path):
    if not bin_path("pga-hts", "bcftools").exists():
        pytest.skip("pga-hts no disponible")
    sources = _make_sources(tmp_path / "src")
    calib = _load_script()
    report = calib.main([
        "--workdir", str(tmp_path / "work"),
        "--sources", str(sources),
        "--ref-template", str(tmp_path / "src" / "ucsc" / "chr{chrom}.fa.gz"),
        "--onekg-base", str(tmp_path / "src" / "onekg"),
        "--ks", "50,100,200", "--reps", "40",
    ])
    work = tmp_path / "work"
    # contigs renombrados, SM reetiquetado, panel sin sitios raros ni sin cobertura
    assert (work / "ref" / "grch37_subset.fa").read_text().startswith(">21\n")
    assert pysam.AlignmentFile(str(work / "bams" / "HG002.bam")).header["RG"][0]["SM"] == "HG002"
    sites = [l.split("\t") for l in (work / "panel" / "sites.tsv").read_text().splitlines()]
    assert sites and all(c in ("21", "22") for c, *_ in sites)
    cal = report["calibration"]["scenarios"]
    true_sc = cal["trío verdadero (HG003, HG004 -> HG002)"]
    false_sc = cal["padre falso (HG005 declarado padre de HG002)"]
    assert true_sc["all_sites_status"] == "passed"
    assert false_sc["all_sites_status"] == "failed"
    row200 = {r["k"]: r for r in false_sc["by_panel_size"]}[200]
    assert row200["correct_rate"] > 0.95
    assert (work / "report.md").read_text().startswith("# Calibración")
