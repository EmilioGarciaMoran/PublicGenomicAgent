"""pga run-family de punta a punta: familia consanguínea con miembros sin BAM."""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

from publicgenomicagent.env.micromamba import bin_path

sys.path.insert(0, str(Path(__file__).parents[1] / "fixtures"))
pysam = pytest.importorskip("pysam")
pytestmark = pytest.mark.integration

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


def _gts(vcf: Path):
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = c[9:]
                continue
            fmt = c[8].split(":")
            yield c, {s: dict(zip(fmt, v.split(":"))) for s, v in zip(samples, c[9:])}


def test_run_family_consanguineous(tmp_path, monkeypatch):
    if not bin_path("pga-hts", "bcftools").exists():
        pytest.skip("pga-hts no disponible")
    from make_family_bams import make_family
    from publicgenomicagent.cli.main import app

    fam = make_family(tmp_path / "fam", PARENTS, n_sites=300, seed=21)
    (tmp_path / "fam.ped").write_text(PED)
    monkeypatch.chdir(tmp_path)
    args = ["run-family", "--ped", str(tmp_path / "fam.ped"),
            "--reference", str(fam["ref"]), "--region", "chr1:400-60500",
            "--case-id", "consang", "--relatedness-sites", str(fam["sites"])]
    for s in SEQ:
        args += ["--bam", f"{s}={fam['bams'][s]}"]
    ambiguous = CliRunner().invoke(app, args)          # dos afectados: no se adivina
    assert ambiguous.exit_code == 2 and "--proband" in ambiguous.output
    r = CliRunner().invoke(app, args + ["--proband", "P1"], env={"COLUMNS": "200"})
    assert r.exit_code == 0, r.output

    st = json.loads((tmp_path / "results" / "consang.session.json").read_text())
    tools = [c["tool_name"] for c in st["tool_calls"]]
    assert tools.count("qc_bam") == 5 and "verify_relatedness" in tools
    assert "joint_call" in tools and "mendelian_filter" in tools
    rel = next(c for c in st["tool_calls"] if c["tool_name"] == "verify_relatedness")["output"]
    assert rel["status"] in ("passed", "warning"), rel["issues"]   # consanguinidad sin falsos fallos

    mf = next(c for c in st["tool_calls"] if c["tool_name"] == "mendelian_filter")["output"]
    rep = json.loads(Path(mf["report_json"]).read_text())
    assert rep["family_roles"] == {"affected": ["P1", "P2"], "carriers": ["FA", "MO"],
                                   "carriers_any": [], "unaffected": ["S1"]}
    # recuento independiente sobre el VCF conjunto
    joint = next(c for c in st["tool_calls"] if c["tool_name"] == "joint_call")["output"]["output_vcf"]
    hom = lambda g: g.get("GT") in ("1/1", "1|1")
    het = lambda g: g.get("GT") in ("0/1", "1/0", "0|1", "1|0")
    dp_ok = lambda g: g.get("DP", "0").isdigit() and int(g["DP"]) >= 10
    expected = sum(
        1 for c, g in _gts(Path(joint))
        if float(c[5]) >= 20 and all(hom(g[a]) and dp_ok(g[a]) for a in ("P1", "P2"))
        and het(g["FA"]) and het(g["MO"]) and not hom(g["S1"]))
    assert rep["counts"]["auto_rec_hom"] == expected
    assert "Segregación" in r.output and "sin BAM ['GF', 'GM', 'A', 'B', 'SA', 'SB']" in r.output


def test_run_family_compound_het_consistent_with_joint_vcf(tmp_path, monkeypatch):
    """Cada pareja de heterocigoto compuesto informada cumple las reglas sobre el VCF real."""
    if not bin_path("pga-hts", "bcftools").exists():
        pytest.skip("pga-hts no disponible")
    from make_family_bams import make_family
    from publicgenomicagent.cli.main import app

    fam = make_family(tmp_path / "fam", PARENTS, n_sites=300, seed=21)
    (tmp_path / "fam.ped").write_text(PED)
    monkeypatch.chdir(tmp_path)
    args = ["run-family", "--ped", str(tmp_path / "fam.ped"), "--reference", str(fam["ref"]),
            "--region", "chr1:400-60500", "--label", "GENSIM", "--case-id", "ch", "--proband", "P1"]
    for s in SEQ:
        args += ["--bam", f"{s}={fam['bams'][s]}"]
    r = CliRunner().invoke(app, args, env={"COLUMNS": "200"})
    assert r.exit_code == 0, r.output
    st = json.loads((tmp_path / "results" / "ch.session.json").read_text())
    mf = next(c for c in st["tool_calls"] if c["tool_name"] == "mendelian_filter")
    assert mf["input"]["regions"][0]["label"] == "GENSIM"
    rep = json.loads(Path(mf["output"]["report_json"]).read_text())
    assert "comp_het" in rep["counts"] and "comp_het_unphased" in rep["counts"]
    joint = next(c for c in st["tool_calls"] if c["tool_name"] == "joint_call")["output"]["output_vcf"]
    g = {f"{c[0]}\t{c[1]}": gt for c, gt in _gts(Path(joint))}
    het = lambda x: x.get("GT") in ("0/1", "1/0", "0|1", "1|0")
    carrier = lambda x: x.get("GT") in ("0/1", "1/0", "0|1", "1|0", "1/1", "1|1")
    ref = lambda x: x.get("GT") in ("0/0", "0|0")
    assert rep["comp_het_pairs"], "la familia simulada debería tener parejas candidatas"
    for p in rep["comp_het_pairs"]:
        a, b = (g["\t".join(p[k].split("\t")[:2])] for k in ("v1", "v2"))
        assert all(het(a[s]) and het(b[s]) for s in ("P1", "P2"))
        assert not any(carrier(a[s]) and carrier(b[s]) for s in ("FA", "MO", "S1"))
        informative = [s for s in ("FA", "MO", "S1")
                       if (carrier(a[s]) and ref(b[s])) or (carrier(b[s]) and ref(a[s]))]
        assert p["phase"] == ("trans" if informative else "unknown")
        assert sorted(p["phase_evidence"]) == sorted(informative)
    assert rep["counts"]["comp_het"] == sum(p["phase"] == "trans" for p in rep["comp_het_pairs"])


def test_run_family_rare_variant_filter(tmp_path, monkeypatch):
    """--max-af: con gnomAD simulado, ninguna variante común llega a los modelos."""
    if not bin_path("pga-hts", "bcftools").exists():
        pytest.skip("pga-hts no disponible")
    import publicgenomicagent.tools.frequency as freq
    from make_family_bams import make_family
    from publicgenomicagent.cli.main import app

    fam = make_family(tmp_path / "fam", PARENTS, n_sites=300, seed=21)
    (tmp_path / "fam.ped").write_text(PED)
    sites = [l.split("\t") for l in Path(fam["sites"]).read_text().splitlines()]
    common = {int(p) for i, (_, p, _, _) in enumerate(sites) if i % 10}   # 90% comunes
    mid_only = {int(sites[0][1])}                                         # común solo en MID

    def fake_fetch(chrom, start, stop, refresh=False):
        out = []
        for c, p, r, a in sites:
            p = int(p)
            if not start <= p <= stop:
                continue
            if p in mid_only:
                ex = {"ac": 1, "an": 10000, "populations": [{"id": "mid", "ac": 30, "an": 1000}]}
            elif p in common:
                ex = {"ac": 3000, "an": 10000, "populations": [{"id": "mid", "ac": 300, "an": 1000}]}
            else:
                ex = {"ac": 2, "an": 10000, "populations": [{"id": "mid", "ac": 0, "an": 1000}]}
            out.append({"variant_id": f"1-{p}-{r}-{a}", "pos": p, "ref": r, "alt": a,
                        "exome": ex, "genome": None})
        return out

    monkeypatch.setattr(freq, "fetch_region_af", fake_fetch)
    monkeypatch.chdir(tmp_path)
    args = ["run-family", "--ped", str(tmp_path / "fam.ped"), "--reference", str(fam["ref"]),
            "--region", "chr1:400-60500", "--label", "GENSIM", "--case-id", "rare",
            "--proband", "P1", "--max-af", "0.01"]
    for s in SEQ:
        args += ["--bam", f"{s}={fam['bams'][s]}"]
    r = CliRunner().invoke(app, args, env={"COLUMNS": "200"})
    assert r.exit_code == 0, r.output
    st = json.loads((tmp_path / "results" / "rare.session.json").read_text())
    tools = [c["tool_name"] for c in st["tool_calls"]]
    assert tools.index("joint_call") < tools.index("annotate_frequency") < tools.index("mendelian_filter")
    assert st["params"]["max_af"] == 0.01
    mf = next(c for c in st["tool_calls"] if c["tool_name"] == "mendelian_filter")
    assert mf["input"]["max_af"] == 0.01 and mf["input"]["trio_vcf"].endswith(".freq.vcf.gz")
    rep = json.loads(Path(mf["output"]["report_json"]).read_text())
    rf = rep["rarity_filter"]
    assert rf["removed_common"] > 0 and rf["removed_common_mid_only"] == 1
    reported = {int(p[k].split("\t")[1]) for p in rep["comp_het_pairs"] for k in ("v1", "v2")}
    assert not reported & (common | mid_only)
    assert rep["counts"]["comp_het"] + rep["counts"]["comp_het_unphased"] < 224
