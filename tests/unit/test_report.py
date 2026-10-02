"""Tests for the HTML report generator."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from publicgenomicagent.agent.report import render_session_report


def _minimal_session() -> dict:
    return {
        "case": {
            "case_id": "TEST-001",
            "source": "clinical",
            "proband": "proband",
            "pedigree": [
                {"sample": "father", "sex": "M", "affected": False},
                {"sample": "mother", "sex": "F", "affected": False},
                {"sample": "proband", "sex": "U", "affected": True,
                 "father": "father", "mother": "mother"},
            ],
            "candidate_rois": [
                {"chrom": "chr1", "start": 1000, "end": 1200, "label": "DEMO1"},
            ],
        },
        "tool_calls": [
            {"tool_name": "qc_bam", "ok": True, "duration_ms": 5,
             "input": {}, "output": {}, "error": None, "timestamp": "2026-01-01T00:00:00Z"},
        ],
        "vcfs": {},
        "outputs": {"mendelian_filter": {"counts": {"auto_rec_hom": 1}}},
        "bams": {}, "references": {}, "artifacts": {}, "notes": [],
    }


def test_render_report_creates_html(tmp_path: Path):
    session = tmp_path / "session.json"
    session.write_text(json.dumps(_minimal_session()))
    out = tmp_path / "report.html"

    result = render_session_report(session, out)
    assert result == out
    assert out.exists()

    html = out.read_text()
    assert "<!DOCTYPE html>" in html
    assert "TEST-001" in html
    assert "proband" in html
    assert "qc_bam" in html


def test_render_report_handles_empty_case(tmp_path: Path):
    session = tmp_path / "session.json"
    session.write_text(json.dumps({"case": {"case_id": "EMPTY"}}))
    out = tmp_path / "report.html"

    render_session_report(session, out)
    html = out.read_text()
    assert "EMPTY" in html
    # No debe contener rutas absolutas al LLM (por privacidad)
    # (el report sí puede contener paths, es local)
    assert "PublicGenomicAgent" in html


def test_render_report_missing_session_raises(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        render_session_report(tmp_path / "nope.json", tmp_path / "out.html")

# ---------------------------------------------------------------------------
# Sección de variantes (lectura de VCFs)
# ---------------------------------------------------------------------------

def test_read_vcf_variants_basic(tmp_path: Path):
    """_read_vcf_variants lee un VCF.gz mínimo."""
    import gzip
    from publicgenomicagent.agent.report import _read_vcf_variants

    vcf_path = tmp_path / "test.vcf.gz"
    with gzip.open(vcf_path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr1,length=10000>\n")
        f.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar significance">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfather\tmother\tproband\n")
        f.write("chr1\t1100\t.\tA\tG\t540.8\tPASS\tCLNSIG=Pathogenic\tGT:DP\t0/1:63\t0/1:63\t1/1:63\n")

    variants = _read_vcf_variants(vcf_path)
    assert len(variants) == 1
    v = variants[0]
    assert v["chrom"] == "chr1"
    assert v["pos"] == "1100"
    assert v["ref"] == "A"
    assert v["alt"] == "G"
    assert v["genotypes"] == {"father": "0/1", "mother": "0/1", "proband": "1/1"}
    assert v["CLNSIG"] == "Pathogenic"


def test_render_variants_table_with_annotations():
    """_render_variants_table maneja variantes con CLNSIG."""
    from publicgenomicagent.agent.report import _render_variants_table

    variants = [
        {
            "chrom": "chr1", "pos": "1100", "ref": "A", "alt": "G",
            "qual": "540.8", "filter": "PASS",
            "genotypes": {"father": "0/1", "mother": "0/1", "proband": "1/1"},
            "CLNSIG": "Pathogenic", "CLNDN": "Test disease",
        },
    ]
    html = _render_variants_table(variants, ["father", "mother", "proband"])
    assert "chr1:1100" in html
    assert "A&gt;G" in html
    assert "0/1" in html
    assert "1/1" in html
    assert "Pathogenic" in html
    assert "Test disease" in html


def test_render_variants_table_empty():
    """_render_variants_table maneja lista vacía."""
    from publicgenomicagent.agent.report import _render_variants_table

    html = _render_variants_table([], ["father"])
    assert "No variants" in html

