from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from publicgenomicagent.tools.base import LocalPangenomeInput
from publicgenomicagent.tools.local_pangenome import local_pangenome


FIXTURES = Path(__file__).parents[1] / "fixtures"
MAKE_COHORT = FIXTURES / "make_cohort_fixture.py"
REF_CACHE = Path("~/.pga/cache/reference").expanduser()
BCFTOOLS = str(Path("~/.pga/envs/pga-hts/bin/bcftools").expanduser())
SAMTOOLS = str(Path("~/.pga/envs/pga-hts/bin/samtools").expanduser())

ROI = "chr2_roi:1-25000"
CONTIG = "chr2_roi"


def _env_ready() -> bool:
    return Path(BCFTOOLS).exists() and Path(SAMTOOLS).exists()


def _build_reference(out_fasta: Path) -> None:
    """Extrae la secuencia del ROI y la renombra a chr2_roi."""
    src = REF_CACHE / "hg38_chr2_110000001_110025000.fa"
    if not src.exists():
        pytest.skip(f"referencia cache no disponible: {src}")
    lines = src.read_text().splitlines()
    seq = "".join(l for l in lines if not l.startswith(">"))
    with out_fasta.open("w") as f:
        f.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
    subprocess.run([SAMTOOLS, "faidx", str(out_fasta)], check=True)


@pytest.fixture
def workdir(tmp_path: Path) -> dict:
    cohort_dir = tmp_path / "cohort"
    cohort_dir.mkdir()
    subprocess.run(
        [sys.executable, str(MAKE_COHORT), str(cohort_dir)],
        check=True, capture_output=True,
    )
    ref = cohort_dir / "nphp1_ref.fa"
    _build_reference(ref)
    return {
        "tmp": tmp_path,
        "cohort": cohort_dir / "cohort.vcf.gz",
        "reference": ref,
        "outdir": tmp_path / "pangenome",
    }


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_local_pangenome_produces_three_outputs(workdir: dict) -> None:
    result = local_pangenome(LocalPangenomeInput(
        cohort_vcf=workdir["cohort"],
        reference_fasta=workdir["reference"],
        region=ROI,
        output_dir=workdir["outdir"],
        min_af=0.05,
        af_info_field="AF_MID",
        contig_name=CONTIG,
    ))

    assert result.ok
    assert result.enriched_fasta.exists()
    assert result.enriched_fai.exists()
    assert result.diff_vcf.exists()
    assert result.diff_tbi.exists()
    assert result.report_json.exists()


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_local_pangenome_enriched_differs_from_ref(workdir: dict) -> None:
    result = local_pangenome(LocalPangenomeInput(
        cohort_vcf=workdir["cohort"],
        reference_fasta=workdir["reference"],
        region=ROI,
        output_dir=workdir["outdir"],
        min_af=0.05,
        af_info_field="AF_MID",
        contig_name=CONTIG,
    ))

    def seq(p: Path) -> str:
        return "".join(
            l.strip() for l in p.read_text().splitlines() if not l.startswith(">")
        )

    ref = seq(workdir["reference"])
    enr = seq(result.enriched_fasta)
    assert ref != enr
    # Longitud muy similar (indels pequeños), pero no idéntica
    assert abs(len(enr) - len(ref)) < 200


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_local_pangenome_af_threshold_filters(workdir: dict) -> None:
    """Umbral alto => menos variantes aplicadas."""
    low = local_pangenome(LocalPangenomeInput(
        cohort_vcf=workdir["cohort"],
        reference_fasta=workdir["reference"],
        region=ROI,
        output_dir=workdir["tmp"] / "low",
        min_af=0.01,
        af_info_field="AF_MID",
        contig_name=CONTIG,
    ))
    high = local_pangenome(LocalPangenomeInput(
        cohort_vcf=workdir["cohort"],
        reference_fasta=workdir["reference"],
        region=ROI,
        output_dir=workdir["tmp"] / "high",
        min_af=0.20,
        af_info_field="AF_MID",
        contig_name=CONTIG,
    ))
    assert low.counts["after_af_filter"] >= high.counts["after_af_filter"]
    assert high.counts["after_af_filter"] < low.counts["after_af_filter"]


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_local_pangenome_report_has_hashes(workdir: dict) -> None:
    result = local_pangenome(LocalPangenomeInput(
        cohort_vcf=workdir["cohort"],
        reference_fasta=workdir["reference"],
        region=ROI,
        output_dir=workdir["outdir"],
        min_af=0.05,
        af_info_field="AF_MID",
        contig_name=CONTIG,
    ))

    report = json.loads(result.report_json.read_text())
    assert "cohort_vcf_sha256" in report
    assert "enriched_fasta_sha256" in report
    assert report["counts"]["variants_applied"] > 0


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_local_pangenome_ref_mismatch_aborts(workdir: dict) -> None:
    """Si el FASTA declara otro contig y no hay override, debe abortar."""
    with pytest.raises(RuntimeError, match="REF mismatch|FASTA declara"):
        local_pangenome(LocalPangenomeInput(
            cohort_vcf=workdir["cohort"],
            reference_fasta=workdir["reference"],
            region=ROI,
            output_dir=workdir["tmp"] / "mismatch",
            min_af=0.05,
            af_info_field="AF_MID",
            contig_name="chr_no_existe",
        ))
