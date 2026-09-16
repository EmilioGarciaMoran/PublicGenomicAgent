from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import pytest

from publicgenomicagent.env.paths import registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.tools.base import CallVariantsInput
from publicgenomicagent.tools.call_variants import call_variants


FIXTURES = Path(__file__).parents[1] / "fixtures"
MAKE_BAM = FIXTURES / "make_tiny_bam.sh"
MAKE_REF = FIXTURES / "make_tiny_ref.py"


def _env_ready() -> bool:
    try:
        reg = load_registry(registry_path())
        runtime = ToolRuntime(reg)
        runtime.version("samtools")
        runtime.version("bcftools")
        return True
    except Exception:
        return False


def _read_vcf_lines(vcf_gz: Path) -> list[str]:
    with gzip.open(vcf_gz, "rt") as f:
        return [l for l in f if not l.startswith("#")]


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_call_variants_produces_indexed_vcf(tmp_path: Path) -> None:
    runtime = ToolRuntime(load_registry(registry_path()))
    samtools = str(Path("~/.pga/envs/pga-hts/bin/samtools").expanduser())

    bam = tmp_path / "tiny.bam"
    ref = tmp_path / "ref.fa"
    subprocess.run([str(MAKE_BAM), samtools, str(bam), "CHILD"], check=True)
    subprocess.run(["python3", str(MAKE_REF), str(ref)], check=True)

    out = tmp_path / "child.vcf.gz"
    result = call_variants(
        runtime,
        CallVariantsInput(
            bam_path=bam,
            reference_fasta=ref,
            output_vcf=out,
            min_qual=0,
            min_dp=1,
        ),
    )

    assert result.ok
    assert result.output_vcf.exists()
    assert result.output_tbi.exists()
    assert result.variants_total > 0
    assert result.variants_total == result.variants_passing

    lines = _read_vcf_lines(result.output_vcf)
    assert len(lines) == result.variants_total
    # Todas las variantes deben estar en chr1
    for line in lines:
        assert line.split("\t")[0] == "chr1"


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_call_variants_rejects_missing_bam(tmp_path: Path) -> None:
    runtime = ToolRuntime(load_registry(registry_path()))
    ref = tmp_path / "ref.fa"
    subprocess.run(["python3", str(MAKE_REF), str(ref)], check=True)

    with pytest.raises(FileNotFoundError):
        call_variants(
            runtime,
            CallVariantsInput(
                bam_path=tmp_path / "noexiste.bam",
                reference_fasta=ref,
                output_vcf=tmp_path / "out.vcf.gz",
            ),
        )


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_call_variants_region_filter_narrows_variants(tmp_path: Path) -> None:
    """Llamar solo chr1:100-300 debe producir menos variantes que el BAM completo."""
    runtime = ToolRuntime(load_registry(registry_path()))
    samtools = str(Path("~/.pga/envs/pga-hts/bin/samtools").expanduser())

    bam = tmp_path / "tiny.bam"
    ref = tmp_path / "ref.fa"
    subprocess.run([str(MAKE_BAM), samtools, str(bam), "CHILD"], check=True)
    subprocess.run(["python3", str(MAKE_REF), str(ref)], check=True)

    full = call_variants(
        runtime,
        CallVariantsInput(
            bam_path=bam, reference_fasta=ref,
            output_vcf=tmp_path / "full.vcf.gz",
            min_qual=0, min_dp=1,
        ),
    )
    narrowed = call_variants(
        runtime,
        CallVariantsInput(
            bam_path=bam, reference_fasta=ref,
            output_vcf=tmp_path / "roi.vcf.gz",
            region="chr1:100-300",
            min_qual=0, min_dp=1,
        ),
    )

    assert narrowed.variants_total < full.variants_total
    assert narrowed.variants_total > 0
