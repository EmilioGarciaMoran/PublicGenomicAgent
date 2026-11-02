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
