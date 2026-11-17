"""run-family --hpo: LIRICAL en el flujo familiar y cruce fenotipo × genotipo en el informe.

LIRICAL se sustituye por un doble que escribe un ranking fijo (sus datos no se
descargan en CI); el resto del flujo es real.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

from publicgenomicagent.env.micromamba import bin_path

sys.path.insert(0, str(Path(__file__).parents[1] / "fixtures"))
pytest.importorskip("pysam")
pytestmark = pytest.mark.integration
MINI = Path(__file__).parents[1] / "fixtures" / "lirical_mini"

from tests.integration.test_run_family import PARENTS, PED, SEQ  # noqa: E402


def test_run_family_hpo_concordance_in_report(tmp_path, monkeypatch):
    if not bin_path("pga-hts", "bcftools").exists():
        pytest.skip("pga-hts no disponible")
    import publicgenomicagent.tools.registry as reg
    from make_family_bams import make_family
    from publicgenomicagent.cli.main import app
    from publicgenomicagent.tools.base import PhenotypeRankingOutput

    seen = {}

    def fake_lirical(inp):
        seen["hpo"] = list(inp.hpo_ids)
        out = Path(inp.output_dir); out.mkdir(parents=True, exist_ok=True)
        tsv = out / "lirical.tsv"; shutil.copy(MINI / "ranking.tsv", tsv)
        rep = out / "report.json"; rep.write_text("{}")
        return PhenotypeRankingOutput(ok=True, tool="phenotype_ranking", message="doble",
                                      ranking_tsv=tsv, report_json=rep, data_dir=MINI)

    import dataclasses
    monkeypatch.setitem(reg.TOOL_REGISTRY, "phenotype_ranking",
                        dataclasses.replace(reg.TOOL_REGISTRY["phenotype_ranking"], func=fake_lirical))
    fam = make_family(tmp_path / "fam", PARENTS, n_sites=300, seed=21)
    (tmp_path / "fam.ped").write_text(PED)
    monkeypatch.chdir(tmp_path)
    args = ["run-family", "--ped", str(tmp_path / "fam.ped"), "--reference", str(fam["ref"]),
            "--region", "chr1:400-60500", "--label", "NPHP1", "--case-id", "pg", "--proband", "P1",
            "--hpo", "HP:0000090,HP:0000083,HP:0001903"]
    for s in SEQ:
        args += ["--bam", f"{s}={fam['bams'][s]}"]
    bad_args = list(args); bad_args[bad_args.index("--hpo") + 1] = "HP:90,x"
    bad = CliRunner().invoke(app, bad_args)
    assert bad.exit_code == 2 and "no válidos" in bad.output
    r = CliRunner().invoke(app, args, env={"COLUMNS": "200"})
    assert r.exit_code == 0, r.output
    assert seen["hpo"] == ["HP:0000090", "HP:0000083", "HP:0001903"]

    session = tmp_path / "results" / "pg.session.json"
    st = json.loads(session.read_text())
    assert "phenotype_ranking" in [c["tool_name"] for c in st["tool_calls"]]
    mf = next(c for c in st["tool_calls"] if c["tool_name"] == "mendelian_filter")["output"]
    rep = json.loads(Path(mf["report_json"]).read_text())
    assert rep["auto_rec_hom_candidates"] and all(c["gene"] == "NPHP1" for c in rep["auto_rec_hom_candidates"])

    out = tmp_path / "r.html"
    rr = CliRunner().invoke(app, ["report", str(session), "--out", str(out), "--viewer", "none"])
    assert rr.exit_code == 0, rr.output
    html = out.read_text()
    sec = html.split("Phenotype–genotype concordance", 1)[1]
    assert "Nephronophthisis 1" in sec and "OMIM:256100" in sec and "4 of 8525 (tied with 3)" in sec
    assert "not combined into one score" in sec
    assert "NPHP4" in sec          # enfermedad bien situada sin candidata en su gen
