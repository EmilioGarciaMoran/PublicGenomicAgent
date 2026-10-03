"""End-to-end test of the local pangenome pipeline.

Runs the four pangenomic tools in sequence on the tiny reference
used by the DEMO_TRIO fixture:

  local_pangenome -> build_local_graph -> align_to_graph -> call_from_graph

The test verifies that each tool produces its expected artefact
(.fa, .vg, .xg, .gam, .vcf.gz) and that no step crashes. It does
NOT verify the biological content of the variants: that is covered
by the trio pipeline tests (test_demo_trio.py).

Marked as `integration`: requires pga-hts and pga-pangenome.
"""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import pytest

from publicgenomicagent.env.micromamba import bin_path
from publicgenomicagent.env.paths import registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.tools.base import (
    AlignToGraphInput,
    BuildLocalGraphInput,
    CallFromGraphInput,
    LocalPangenomeInput,
)
from publicgenomicagent.tools.local_graph import (
    align_to_graph,
    build_local_graph,
    call_from_graph,
)
from publicgenomicagent.tools.local_pangenome import local_pangenome


pytestmark = pytest.mark.integration


DEMO_DIR = Path("results/demo_trio")


def _envs_ready() -> bool:
    try:
        reg = load_registry(registry_path())
        rt = ToolRuntime(reg)
        rt.version("vg")
        rt.version("samtools")
        rt.version("bcftools")
        return True
    except Exception:  # noqa: BLE001
        return False


@pytest.fixture(scope="module")
def tiny_ref(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Copia la tiny ref del demo o la genera si no existe."""
    out_dir = tmp_path_factory.mktemp("pangenome_e2e")
    ref_src = DEMO_DIR / "ref.fa"
    fai_src = DEMO_DIR / "ref.fa.fai"

    if not ref_src.exists():
        pytest.skip(f"tiny reference not found: {ref_src}")

    ref = out_dir / "ref.fa"
    ref.write_bytes(ref_src.read_bytes())
    (out_dir / "ref.fa.fai").write_bytes(fai_src.read_bytes())
    return ref


@pytest.fixture(scope="module")
def tiny_cohort(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Genera una cohort VCF mínima con una variante en chr1:1100."""
    out_dir = tmp_path_factory.mktemp("pangenome_e2e")
    vcf = out_dir / "cohort.vcf.gz"
    bcftools = bin_path("pga-hts", "bcftools")

    plain = out_dir / "cohort.vcf"
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=10000>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1100\t.\tA\tG\t.\t.\tAF=0.25\n"
    )
    subprocess.run(
        [str(bcftools), "view", "-Oz", "-o", str(vcf), str(plain)],
        check=True, capture_output=True,
    )
    subprocess.run([str(bcftools), "index", "-t", str(vcf)], check=True)
    plain.unlink()
    return vcf


@pytest.mark.skipif(not _envs_ready(), reason="pga-pangenome not available")
def test_pangenome_pipeline_end_to_end(
    tiny_ref: Path,
    tiny_cohort: Path,
    tmp_path: Path,
) -> None:
    reg = load_registry(registry_path())
    rt = ToolRuntime(reg)

    # 1. local_pangenome
    pg_dir = tmp_path / "pangenome"
    pg_inp = LocalPangenomeInput(
        cohort_vcf=tiny_cohort,
        reference_fasta=tiny_ref,
        region="chr1:1-10000",
        output_dir=pg_dir,
        min_af=0.01,
    )
    pg_out = local_pangenome(pg_inp)
    assert pg_out.ok
    assert pg_out.enriched_fasta.exists(), pg_out.enriched_fasta
    assert pg_out.diff_vcf.exists(), pg_out.diff_vcf

    # 2. build_local_graph
    graph_vg = tmp_path / "graph.vg"
    build_inp = BuildLocalGraphInput(
        reference_fasta=pg_out.enriched_fasta,
        cohort_vcf=tiny_cohort,
        output_vg=graph_vg,
    )
    build_out = build_local_graph(rt, build_inp)
    assert build_out.ok
    assert build_out.graph_vg.exists(), build_out.graph_vg
    assert build_out.graph_xg.exists(), build_out.graph_xg

    # 3. align_to_graph: necesitamos FASTQ. Extraemos del BAM del demo.
    #    Usamos el probando (que tiene la variante).
    proband_bam = DEMO_DIR / "proband.bam"
    if not proband_bam.exists():
        pytest.skip(f"proband.bam not found: {proband_bam}")

    fastq = tmp_path / "proband.fastq"
    samtools = bin_path("pga-hts", "samtools")
    subprocess.run(
        [str(samtools), "fastq", "-0", str(fastq), str(proband_bam)],
        check=True, capture_output=True,
    )
    assert fastq.exists()

    gam = tmp_path / "reads.gam"
    align_inp = AlignToGraphInput(
        graph_vg=build_out.graph_vg,
        reads_fastq=fastq,
        output_gam=gam,
    )
    align_out = align_to_graph(rt, align_inp)
    assert align_out.ok
    assert align_out.gam.exists(), align_out.gam

    # 4. call_from_graph: requiere un pack de cobertura.
    #    Lo creamos con vg pack.
    pack = tmp_path / "reads.pack"
    vg_bin = bin_path("pga-pangenome", "vg")
    subprocess.run(
        [str(vg_bin), "pack", "-x", str(build_out.graph_xg),
         "-g", str(align_out.gam), "-o", str(pack)],
        check=True, capture_output=True,
    )
    assert pack.exists()

    call_vcf = tmp_path / "calls.vcf.gz"
    call_inp = CallFromGraphInput(
        graph_vg=build_out.graph_vg,
        graph_xg=build_out.graph_xg,
        pack=pack,
        output_vcf=call_vcf,
    )
    call_out = call_from_graph(rt, call_inp)
    assert call_out.ok
    assert call_out.vcf.exists(), call_out.vcf
    # El VCF debe ser un gzip válido
    with gzip.open(call_out.vcf, "rt") as f:
        header = f.readline()
    assert header.startswith("##fileformat=VCF")
