from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path

import pytest

from publicgenomicagent.tools.base import CompareVCFsInput
from publicgenomicagent.tools.compare_vcfs import compare_vcfs


BCFTOOLS = str(Path("~/.pga/envs/pga-hts/bin/bcftools").expanduser())
FIXTURES = Path(__file__).parents[1] / "fixtures"
MAKE_COHORT = FIXTURES / "make_cohort_fixture.py"


def _env_ready() -> bool:
    try:
        return Path(BCFTOOLS).exists()
    except Exception:
        return False


@pytest.fixture
def cohort_dir(tmp_path: Path) -> Path:
    """Genera el fixture en tmp_path y devuelve el directorio."""
    import sys as _sys
    subprocess.run(
        [_sys.executable, str(MAKE_COHORT), str(tmp_path)],
        check=True, capture_output=True,
    )
    return tmp_path


def _common_snvs(cohort_vcf: Path, out_vcf: Path, min_af: float) -> Path:
    """Filtra cohort.vcf.gz a SNVs+indels con INFO/AF >= min_af."""
    # Paso 1: filtrar por AF (pero esto deja pasar SVs que no tienen AF)
    tmp = out_vcf.with_suffix(".tmp.vcf.gz")
    subprocess.run(
        [BCFTOOLS, "view", "-i", f"INFO/AF>={min_af}",
         str(cohort_vcf), "-Oz", "-o", str(tmp)],
        check=True,
    )
    # Paso 2: excluir SVs simbólicos
    subprocess.run(
        [BCFTOOLS, "view", "-v", "snps,indels",
         str(tmp), "-Oz", "-o", str(out_vcf)],
        check=True,
    )
    subprocess.run(
        [BCFTOOLS, "index", "-t", str(out_vcf)],
        check=True,
    )
    tmp.unlink(missing_ok=True)
    tmp.with_suffix(".tmp.vcf.gz.tbi").unlink(missing_ok=True)
    return out_vcf


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_compare_vcfs_self_is_all_consistent(cohort_dir: Path) -> None:
    """Comparar un VCF consigo mismo: todo consistent, nada recovered/lost."""
    vcf = cohort_dir / "cohort.vcf.gz"
    out = cohort_dir / "delta_self.vcf.gz"
    report = cohort_dir / "delta_self.json"

    result = compare_vcfs(CompareVCFsInput(
        baseline_vcf=vcf,
        candidate_vcf=vcf,
        output_delta_vcf=out,
        output_report=report,
        sample="C2",
    ))

    assert result.counts.get("consistent", 0) > 0
    assert result.counts.get("recovered", 0) == 0
    assert result.counts.get("lost", 0) == 0
    assert result.counts.get("discordant", 0) == 0


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_compare_vcfs_subset_is_recovered(cohort_dir: Path) -> None:
    """Baseline subconjunto -> todo recovered o consistent, nada lost."""
    cohort = cohort_dir / "cohort.vcf.gz"
    baseline = cohort_dir / "baseline_common.vcf.gz"
    candidate = cohort_dir / "candidate_all.vcf.gz"

    # Construir baseline (AF>=0.05, SNVs+indels)
    _common_snvs(cohort, baseline, min_af=0.05)

    # Construir candidate (todas SNVs+indels, sin SVs)
    subprocess.run(
        [BCFTOOLS, "view", "-v", "snps,indels",
         str(cohort), "-Oz", "-o", str(candidate)],
        check=True,
    )
    subprocess.run([BCFTOOLS, "index", "-t", str(candidate)], check=True)

    result = compare_vcfs(CompareVCFsInput(
        baseline_vcf=baseline,
        candidate_vcf=candidate,
        output_delta_vcf=cohort_dir / "delta_af.vcf.gz",
        output_report=cohort_dir / "delta_af.json",
        sample="C2",
    ))

    # baseline ⊂ candidate => lost = 0
    assert result.counts.get("lost", 0) == 0
    assert result.counts.get("recovered", 0) > 0
    assert result.counts.get("consistent", 0) > 0

    # Reporte escrito correctamente
    report = json.loads((cohort_dir / "delta_af.json").read_text())
    assert report["baseline_total"] == result.counts.get("consistent", 0)
    assert report["candidate_total"] == (
        result.counts.get("recovered", 0)
        + result.counts.get("consistent", 0)
    )


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_compare_vcfs_delta_vcf_has_status_info(cohort_dir: Path) -> None:
    """El delta.vcf.gz debe llevar INFO/PGA_STATUS en cada variante."""
    cohort = cohort_dir / "cohort.vcf.gz"
    out = cohort_dir / "delta_self.vcf.gz"

    compare_vcfs(CompareVCFsInput(
        baseline_vcf=cohort,
        candidate_vcf=cohort,
        output_delta_vcf=out,
        sample="C2",
    ))

    with gzip.open(out, "rt") as f:
        data_lines = [l for l in f if not l.startswith("#")]
    assert len(data_lines) > 0
    for line in data_lines:
        fields = line.split("\t")
        assert "PGA_STATUS=" in fields[7]
        assert "PGA_BASELINE_GT=" in fields[7]
        assert "PGA_CANDIDATE_GT=" in fields[7]


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-hts creado")
def test_compare_vcfs_missing_file_raises(cohort_dir: Path) -> None:
    with pytest.raises(FileNotFoundError):
        compare_vcfs(CompareVCFsInput(
            baseline_vcf=cohort_dir / "noexiste.vcf.gz",
            candidate_vcf=cohort_dir / "cohort.vcf.gz",
            output_delta_vcf=cohort_dir / "out.vcf.gz",
        ))
