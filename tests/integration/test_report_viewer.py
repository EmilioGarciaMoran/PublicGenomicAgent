"""E3: informe con verificación primero y visor JBrowse embebido que abre sin red."""
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

from tests.integration.test_run_family import PARENTS, PED, SEQ  # noqa: E402


@pytest.fixture(scope="module")
def session(tmp_path_factory):
    if not bin_path("pga-hts", "bcftools").exists():
        pytest.skip("pga-hts no disponible")
    from make_family_bams import make_family
    from publicgenomicagent.cli.main import app
    tmp = tmp_path_factory.mktemp("rv")
    fam = make_family(tmp / "fam", PARENTS, n_sites=300, seed=21)
    (tmp / "fam.ped").write_text(PED)
    import os
    cwd = os.getcwd(); os.chdir(tmp)
    try:
        args = ["run-family", "--ped", str(tmp / "fam.ped"), "--reference", str(fam["ref"]),
                "--region", "chr1:400-60500", "--case-id", "rv", "--proband", "P1",
                "--relatedness-sites", str(fam["sites"])]
        for s in SEQ:
            args += ["--bam", f"{s}={fam['bams'][s]}"]
        r = CliRunner().invoke(app, args)
        assert r.exit_code == 0, r.output
    finally:
        os.chdir(cwd)
    return tmp / "results" / "rv.session.json"


def _report(session, out, *extra):
    from publicgenomicagent.cli.main import app
    r = CliRunner().invoke(app, ["report", str(session), "--out", str(out), *extra])
    assert r.exit_code == 0, r.output
    return out.read_text()


def test_report_verification_before_analysis_and_viewer(session, tmp_path):
    html = _report(session, tmp_path / "r.html")
    i_ver, i_an, i_view = (html.index(t) for t in
                           ("1. Data verification", "2. Analysis", "3. Embedded viewer"))
    assert i_ver < i_an < i_view
    assert "Family verification (local IBS)" in html and "Per-sample QC" in html
    assert "Privacy:" in html and "pga-jbrowse-data" in html
    data = json.loads(html.split('id="pga-jbrowse-data" type="application/json">')[1]
                      .split("</script>")[0].replace("<\\/", "</"))
    assert data["show"] == ["family_vcf"] + [f"bam_{s}" for s in SEQ]


def test_report_no_reads_and_none(session, tmp_path):
    html = _report(session, tmp_path / "n.html", "--no-reads")
    assert "Privacy:" not in html and '"bam_' not in html and "family_vcf" in html
    html = _report(session, tmp_path / "x.html", "--viewer", "none")
    assert "pga-jbrowse-data" not in html and "1. Data verification" in html


def test_viewer_renders_offline_in_browser(session, tmp_path):
    pw = pytest.importorskip("playwright.sync_api")
    out = tmp_path / "r.html"
    _report(session, out)
    external, errors = [], []
    with pw.sync_playwright() as p:
        b = p.chromium.launch()
        pg = b.new_page()
        pg.route("**/*", lambda route: (external.append(route.request.url), route.abort())
                 if not route.request.url.startswith(("file:", "data:", "blob:")) else route.continue_())
        pg.on("pageerror", lambda e: errors.append(str(e)))
        pg.goto(out.as_uri())
        pg.wait_for_selector("#pga-jbrowse canvas, #pga-jbrowse svg", timeout=30_000)
        pg.wait_for_timeout(4000)
        text = pg.inner_text("#pga-jbrowse")
        pg.screenshot(path=str(tmp_path / "viewer.png"), full_page=True)
        b.close()
    assert external == [] and errors == []
    assert "Family genotypes" in text and "P1 (reads)" in text


def test_report_redact_paths(session, tmp_path):
    import re
    root = str(session.parents[1])
    html = _report(session, tmp_path / "s.html", "--redact-paths")
    outside_scripts = re.sub(r"<script\b.*?</script>", "", html, flags=re.S)
    assert root in _report(session, tmp_path / "u.html", "--viewer", "none")
    assert root not in outside_scripts and "vcf.gz" in outside_scripts
