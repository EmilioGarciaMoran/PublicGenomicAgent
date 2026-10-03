"""Tests for the IGV integration module."""
from __future__ import annotations

import json
from pathlib import Path

from publicgenomicagent.agent.igv import (
    build_control_url,
    build_session_from_state,
    build_session_xml,
)


def test_build_control_url_goto():
    url = build_control_url(locus="chr1:1000-1200")
    assert url == "http://localhost:60151/goto?locus=chr1:1000-1200"


def test_build_control_url_load(tmp_path: Path):
    f = tmp_path / "test.bam"
    f.touch()
    url = build_control_url(load_path=f)
    assert "http://localhost:60151/load?file=" in url
    assert "test.bam" in url


def test_build_control_url_custom_port():
    url = build_control_url(locus="chr1:1-100", port=12345)
    assert url.startswith("http://localhost:12345/")


def test_build_session_xml_basic(tmp_path: Path):
    bam = tmp_path / "proband.bam"
    bam.touch()
    vcf = tmp_path / "variants.vcf.gz"
    vcf.touch()

    xml = build_session_xml(
        reference_path=tmp_path / "ref.fa",
        tracks=[("proband.bam", bam), ("variants.vcf.gz", vcf)],
        locus="chr1:1000-1200",
    )

    assert '<?xml version="1.0"' in xml
    assert '<Session genome="hg38" locus="chr1:1000-1200">' in xml
    assert "proband.bam" in xml
    assert "variants.vcf.gz" in xml
    assert "file://" in xml


def test_build_session_from_state(tmp_path: Path):
    # Crear ficheros dummy que existan
    bam = tmp_path / "proband.bam"
    bam.write_bytes(b"")
    vcf = tmp_path / "variants.vcf.gz"
    vcf.write_bytes(b"")
    ref = tmp_path / "ref.fa"
    ref.write_bytes(b"")

    session = {
        "case": {
            "case_id": "TEST",
            "candidate_rois": [{"chrom": "chr1", "start": 1000, "end": 1200}],
        },
        "bams": {"proband": str(bam)},
        "vcfs": {"called": str(vcf)},
        "references": {"hg38": str(ref)},
    }
    session_file = tmp_path / "session.json"
    session_file.write_text(json.dumps(session))

    info = build_session_from_state(session_file)
    assert info["locus"] == "chr1:1000-1200"
    assert "goto?locus=chr1:1000-1200" in info["goto_url"]
    assert len(info["tracks"]) == 2
    assert "proband.bam" in info["xml"]

