from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from publicgenomicagent.env.paths import registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.tools.base import (
    AlignToGraphInput,
    BuildLocalGraphInput,
    CallFromGraphInput,
)
from publicgenomicagent.tools.local_graph import (
    align_to_graph,
    build_local_graph,
    call_from_graph,
)


FIXTURES = Path(__file__).parents[1] / "fixtures"
MAKE_COHORT = FIXTURES / "make_cohort_fixture.py"
MAKE_BAM = FIXTURES / "make_patient_bam.py"
SAMTOOLS = str(Path("~/.pga/envs/pga-hts/bin/samtools").expanduser())
REF_CACHE = Path("~/.pga/cache/reference").expanduser()
CONTIG = "chr2_roi"


def _vg_ready() -> bool:
    try:
        reg = load_registry(registry_path())
        runtime = ToolRuntime(reg)
        runtime.version("vg")
        return True
    except Exception:
        return False


@pytest.fixture
def workdir(tmp_path: Path) -> dict:
    """Genera cohort, referencia y BAM de C2 en tmp_path."""
    cohort_dir = tmp_path / "cohort"
    cohort_dir.mkdir()
    subprocess.run(
        [sys.executable, str(MAKE_COHORT), str(cohort_dir)],
        check=True, capture_output=True,
    )

    # Referencia del ROI con contig chr2_roi
    src = REF_CACHE / "hg38_chr2_110000001_110025000.fa"
    if not src.exists():
        pytest.skip(f"referencia no cacheada: {src}")
    ref = cohort_dir / "nphp1_ref.fa"
    lines = src.read_text().splitlines()
    seq = "".join(l for l in lines if not l.startswith(">"))
    with ref.open("w") as f:
        f.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
    subprocess.run([SAMTOOLS, "faidx", str(ref)], check=True)

    # BAM de C2
    bam = cohort_dir / "c2.bam"
    subprocess.run(
        [sys.executable, str(MAKE_BAM),
         str(cohort_dir / "cohort.vcf.gz"), str(ref),
         "C2", str(bam), "30"],
        check=True, capture_output=True,
    )

    # FASTQ
    fq = cohort_dir / "c2.fastq"
    with fq.open("wb") as out:
        subprocess.run([SAMTOOLS, "fastq", str(bam)], check=True, stdout=out, stderr=subprocess.DEVNULL)

    return {
        "tmp": tmp_path,
        "cohort": cohort_dir / "cohort.vcf.gz",
        "reference": ref,
        "bam": bam,
        "fastq": fq,
    }


@pytest.mark.integration
@pytest.mark.skipif(not _vg_ready(), reason="requiere pga-pangenome creado")
def test_build_local_graph_produces_vg_and_xg(workdir: dict) -> None:
    runtime = ToolRuntime(load_registry(registry_path()))
    out_vg = workdir["tmp"] / "roi.vg"
    result = build_local_graph(runtime, BuildLocalGraphInput(
        reference_fasta=workdir["reference"],
        cohort_vcf=workdir["cohort"],
        output_vg=out_vg,
    ))
    assert result.ok
    assert result.graph_vg.exists()
    assert result.graph_xg.exists()
    assert result.nodes > 0
    assert result.edges > 0


@pytest.mark.integration
@pytest.mark.skipif(not _vg_ready(), reason="requiere pga-pangenome creado")
def test_align_to_graph_produces_gam(workdir: dict) -> None:
    runtime = ToolRuntime(load_registry(registry_path()))
    out_vg = workdir["tmp"] / "roi.vg"
    build_local_graph(runtime, BuildLocalGraphInput(
        reference_fasta=workdir["reference"],
        cohort_vcf=workdir["cohort"],
        output_vg=out_vg,
    ))
    gam = workdir["tmp"] / "c2.gam"
    result = align_to_graph(runtime, AlignToGraphInput(
        graph_vg=out_vg,
        reads_fastq=workdir["fastq"],
        output_gam=gam,
    ))
    assert result.ok
    assert result.gam.exists()
    assert result.reads_total > 0
    assert result.reads_aligned == result.reads_total


@pytest.mark.integration
@pytest.mark.skipif(not _vg_ready(), reason="requiere pga-pangenome creado")
def test_call_from_graph_produces_vcf(workdir: dict) -> None:
    runtime = ToolRuntime(load_registry(registry_path()))
    out_vg = workdir["tmp"] / "roi.vg"
    build_local_graph(runtime, BuildLocalGraphInput(
        reference_fasta=workdir["reference"],
        cohort_vcf=workdir["cohort"],
        output_vg=out_vg,
    ))
    gam = workdir["tmp"] / "c2.gam"
    pack = workdir["tmp"] / "c2.pack"
    align_to_graph(runtime, AlignToGraphInput(
        graph_vg=out_vg,
        reads_fastq=workdir["fastq"],
        output_gam=gam,
        output_pack=pack,
    ))
    vcf = workdir["tmp"] / "c2.vcf.gz"
    result = call_from_graph(runtime, CallFromGraphInput(
        graph_vg=out_vg,
        graph_xg=out_vg.with_suffix(".xg"),
        pack=pack,
        output_vcf=vcf,
    ))
    assert result.ok
    assert result.vcf.exists()
    assert result.vcf_tbi.exists()
    assert result.variants_total > 0


@pytest.mark.integration
@pytest.mark.skipif(not _vg_ready(), reason="requiere pga-pangenome creado")
def test_call_from_graph_detects_target_variant(workdir: dict) -> None:
    """La variante 2-110008979-A-T debe aparecer en el VCF del grafo."""
    runtime = ToolRuntime(load_registry(registry_path()))
    out_vg = workdir["tmp"] / "roi.vg"
    build_local_graph(runtime, BuildLocalGraphInput(
        reference_fasta=workdir["reference"],
        cohort_vcf=workdir["cohort"],
        output_vg=out_vg,
    ))
    gam = workdir["tmp"] / "c2.gam"
    pack = workdir["tmp"] / "c2.pack"
    align_to_graph(runtime, AlignToGraphInput(
        graph_vg=out_vg,
        reads_fastq=workdir["fastq"],
        output_gam=gam,
        output_pack=pack,
    ))
    vcf = workdir["tmp"] / "c2.vcf.gz"
    call_from_graph(runtime, CallFromGraphInput(
        graph_vg=out_vg,
        graph_xg=out_vg.with_suffix(".xg"),
        pack=pack,
        output_vcf=vcf,
    ))

    # Leer la variante objetivo
    import gzip
    found = False
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[1] == "8979" and fields[3] == "A" and fields[4] == "T":
                found = True
                # Genotipo 0/1
                assert fields[9].startswith("0/1")
                break
    assert found, "variante objetivo 2-110008979-A-T no encontrada"
