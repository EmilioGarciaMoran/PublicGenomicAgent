"""Tabla del clínico -> revisión (dibujo) -> aprobación -> run-family --family."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

from publicgenomicagent.env.micromamba import bin_path

sys.path.insert(0, str(Path(__file__).parents[1] / "fixtures"))
pytest.importorskip("pysam")
pytestmark = pytest.mark.integration

from tests.integration.test_run_family import PARENTS  # noqa: E402

TABLE = """\
# padres primos hermanos; los abuelos y bisabuelos los añade PGA
parents_relationship\tfirst cousins
id\trelation\tsex\tstatus\tbam
P1\tproband\tmale\taffected\tfam/P1.bam
FA\tfather\tmale\tunaffected\tfam/FA.bam
MO\tmother\tfemale\tunaffected\tfam/MO.bam
P2\tsister\tfemale\taffected\tfam/P2.bam
S1\tbrother\tmale\tunaffected\tfam/S1.bam
"""


def test_table_review_approve_run(tmp_path, monkeypatch):
    if not bin_path("pga-hts", "bcftools").exists():
        pytest.skip("pga-hts no disponible")
    from make_family_bams import make_family
    from publicgenomicagent.cli.main import app
    fam = make_family(tmp_path / "fam", PARENTS, n_sites=300, seed=21)
    for s, b in fam["bams"].items():                      # nombres esperados por la tabla
        assert Path(b).name.startswith(s)
    table = tmp_path / "family.tsv"
    table.write_text(TABLE.replace("fam/P1.bam", str(fam["bams"]["P1"])).replace(
        "fam/FA.bam", str(fam["bams"]["FA"])).replace("fam/MO.bam", str(fam["bams"]["MO"])).replace(
        "fam/P2.bam", str(fam["bams"]["P2"])).replace("fam/S1.bam", str(fam["bams"]["S1"])))
    monkeypatch.chdir(tmp_path)
    r = CliRunner()
    run = ["run-family", "--family", str(table), "--reference", str(fam["ref"]),
           "--region", "chr1:400-60500", "--case-id", "tbl", "--relatedness-sites", str(fam["sites"])]
    res = r.invoke(app, run)
    assert res.exit_code == 2 and "sin aprobar" in res.output            # nada corre sin aprobación
    res = r.invoke(app, ["family", "check", str(table)])
    assert res.exit_code == 0, res.output
    assert "F = 1/16" in res.output and (tmp_path / "family.check.html").exists()
    if bin_path("pga-r", "Rscript").exists():
        assert "<svg" in (tmp_path / "family.check.html").read_text()
    res = r.invoke(app, ["family", "approve", str(table), "--by", "Dr Test"])
    assert res.exit_code == 0, res.output
    res = r.invoke(app, run, env={"COLUMNS": "200"})
    assert res.exit_code == 0, res.output
    st = json.loads((tmp_path / "results" / "tbl.session.json").read_text())
    assert st["case"]["proband"] == "P1"
    assert st["params"]["family_approval"]["approved_by"] == "Dr Test"
    rel = next(c for c in st["tool_calls"] if c["tool_name"] == "verify_relatedness")["output"]
    assert rel["status"] in ("passed", "warning"), rel["issues"]       # PED de la tabla = familia real
