"""Test for the annotate_variants tool.

Uses a synthetic annotations VCF with a single variant that
matches the trío DEMO_TRIO (chr1:1100 A>G). Verifies that:

  - bcftools annotate runs without error
  - the output VCF exists and is indexed
  - the copied column (CLNSIG) appears in the output INFO
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from publicgenomicagent.env.micromamba import bin_path
from publicgenomicagent.env.paths import registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.tools.annotate import annotate_variants
from publicgenomicagent.tools.base import AnnotateVariantsInput


pytestmark = pytest.mark.integration


def _vcf_plain_to_bgzf(plain_vcf: Path, bgzf_vcf: Path) -> None:
    """Convierte un VCF de texto plano a VCF bgzipped indexado.

    bcftools index solo acepta ficheros BGZF (no gzip normal). Por eso
    escribimos primero el VCF como texto plano y aquí lo comprimimos
    con `bcftools view -Oz`.
    """
    bcftools = bin_path("pga-hts", "bcftools")
    subprocess.run(
        [str(bcftools), "view", "-Oz", "-o", str(bgzf_vcf), str(plain_vcf)],
        check=True, capture_output=True,
    )
    subprocess.run([str(bcftools), "index", "-t", str(bgzf_vcf)], check=True)
    plain_vcf.unlink()


def _write_annotations_vcf(path: Path) -> None:
    """Escribe un VCF mínimo con una única variante anotada (BGZF)."""
    plain = path.with_suffix(".vcf")
    with plain.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr1,length=10000>\n")
        f.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar significance">\n')
        f.write('##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar disease">\n')
        f.write('##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("chr1\t1100\t.\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNDN=Test+disease;CLNREVSTAT=criteria_provided\n")
    _vcf_plain_to_bgzf(plain, path)


def _write_query_vcf(path: Path) -> None:
    """Escribe un VCF mínimo con la variante a anotar (BGZF)."""
    plain = path.with_suffix(".vcf")
    with plain.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr1,length=10000>\n")
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("chr1\t1100\t.\tA\tG\t100\tPASS\tDP=63\n")
    _vcf_plain_to_bgzf(plain, path)


@pytest.fixture
def runtime() -> ToolRuntime:
    return ToolRuntime(load_registry(registry_path()))


def test_annotate_variants_end_to_end(tmp_path: Path, runtime: ToolRuntime) -> None:
    ann = tmp_path / "annotations.vcf.gz"
    query = tmp_path / "query.vcf.gz"
    out = tmp_path / "annotated.vcf.gz"

    _write_annotations_vcf(ann)
    _write_query_vcf(query)

    inp = AnnotateVariantsInput(
        vcf=query,
        annotations_vcf=ann,
        output_vcf=out,
    )
    result = annotate_variants(runtime, inp)

    assert result.ok
    assert result.variants_total == 1
    assert result.variants_annotated == 1
    assert out.exists()
    assert Path(str(out) + ".tbi").exists()

    # Verificar que el campo CLNSIG está en la salida
    bcftools = bin_path("pga-hts", "bcftools")
    out_text = subprocess.run(
        [str(bcftools), "view", "-H", str(out)],
        capture_output=True, text=True, check=True,
    ).stdout
    assert "CLNSIG=Pathogenic" in out_text
    assert "CLNDN=Test" in out_text


def test_annotate_variants_no_match(tmp_path: Path, runtime: ToolRuntime) -> None:
    """Si la variante no está en las anotaciones, no se anota."""
    ann = tmp_path / "annotations.vcf.gz"
    query = tmp_path / "query.vcf.gz"
    out = tmp_path / "annotated.vcf.gz"

    _write_annotations_vcf(ann)
    # Query con una variante distinta
    plain = query.with_suffix(".vcf")
    with plain.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr1,length=10000>\n")
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("chr1\t5000\t.\tA\tT\t100\tPASS\tDP=63\n")
    _vcf_plain_to_bgzf(plain, query)

    inp = AnnotateVariantsInput(
        vcf=query,
        annotations_vcf=ann,
        output_vcf=out,
    )
    result = annotate_variants(runtime, inp)

    assert result.ok
    assert result.variants_total == 1
    assert result.variants_annotated == 0

