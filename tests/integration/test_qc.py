from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from publicgenomicagent.env.paths import registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.tools.base import QCBamInput, QCLevel
from publicgenomicagent.tools.qc import qc_bam


FIXTURE_SCRIPT = Path(__file__).parents[1] / "fixtures" / "make_tiny_bam.sh"


def _env_ready() -> bool:
    try:
        reg = load_registry(registry_path())
        runtime = ToolRuntime(reg)
        runtime.version("samtools")
        return True
    except Exception:
        return False


@pytest.fixture
def tiny_bam(tmp_path: Path) -> Path:
    """Genera un BAM sintético con sample CHILD."""
    samtools = str(Path("~/.pga/envs/pga-hts/bin/samtools").expanduser())
    out = tmp_path / "tiny.bam"
    subprocess.run(
        [str(FIXTURE_SCRIPT), samtools, str(out), "CHILD"],
        check=True,
    )
    return out


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_qc_structural_passes_on_valid_bam(tiny_bam: Path) -> None:
    runtime = ToolRuntime(load_registry(registry_path()))
    result = qc_bam(runtime, QCBamInput(bam_path=tiny_bam, level=QCLevel.STRUCTURAL))

    assert result.passed
    assert result.issues == []
    assert result.header.sort_order == "coordinate"
    assert result.header.references == ["chr1"]
    assert result.header.samples == ["CHILD"]


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_qc_detects_wrong_expected_sample(tiny_bam: Path) -> None:
    runtime = ToolRuntime(load_registry(registry_path()))
    result = qc_bam(
        runtime,
        QCBamInput(
            bam_path=tiny_bam,
            level=QCLevel.STRUCTURAL,
            expected_sample="WRONG",
        ),
    )
    assert not result.passed
    assert any("WRONG" in i for i in result.issues)


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_qc_counts_produces_flagstat_and_idxstats(tiny_bam: Path) -> None:
    runtime = ToolRuntime(load_registry(registry_path()))
    result = qc_bam(runtime, QCBamInput(bam_path=tiny_bam, level=QCLevel.COUNTS))

    assert result.passed
    assert result.counts is not None
    assert result.counts["idxstats"]["chr1"] == 10
    assert result.counts["flagstat"]["mapped"] == 10


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_qc_missing_bam_reports_issue(tmp_path: Path) -> None:
    runtime = ToolRuntime(load_registry(registry_path()))
    missing = tmp_path / "nope.bam"
    result = qc_bam(runtime, QCBamInput(bam_path=missing))

    assert not result.passed
    assert any("no existe" in i for i in result.issues)
