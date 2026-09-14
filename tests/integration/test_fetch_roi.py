from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from publicgenomicagent.env.paths import registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.tools.base import FetchROIInput
from publicgenomicagent.tools.fetch_roi import fetch_roi


FIXTURE_SCRIPT = Path(__file__).parents[1] / "fixtures" / "make_tiny_bam.sh"


def _env_ready() -> bool:
    try:
        reg = load_registry(registry_path())
        runtime = ToolRuntime(reg)
        runtime.version("samtools")
        return True
    except Exception:
        return False


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_fetch_roi_returns_only_region_reads(tmp_path: Path) -> None:
    reg = load_registry(registry_path())
    runtime = ToolRuntime(reg)

    # 1. Generar BAM sintético usando samtools del entorno pga-hts
    samtools_bin = shutil.which("samtools") or str(
        Path("~/.pga/envs/pga-hts/bin/samtools").expanduser()
    )
    tiny = tmp_path / "tiny.bam"

    import subprocess
    subprocess.run(
        [str(FIXTURE_SCRIPT), samtools_bin, str(tiny)],
        check=True,
    )
    assert tiny.exists()

    # 2. Ejecutar fetch_roi
    roi = tmp_path / "roi.bam"
    out = fetch_roi(
        runtime,
        FetchROIInput(bam_path=tiny, region="chr1:100-300", output_bam=roi),
    )
    assert out.ok
    assert out.output_bam.exists()
    assert out.output_bai.exists()

    # 3. Contar lecturas del sub-BAM
    result = runtime.capture("samtools", ["view", str(out.output_bam)])
    lines = [l for l in result.splitlines() if l.strip()]
    assert len(lines) == 5
    for line in lines:
        pos = int(line.split("\t")[3])
        assert 100 <= pos <= 300
