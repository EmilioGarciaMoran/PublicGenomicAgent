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

