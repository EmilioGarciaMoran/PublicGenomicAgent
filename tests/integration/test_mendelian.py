from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

from publicgenomicagent.env.paths import registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.tools.base import (
    CallVariantsInput,
    IndividualSpec,
    MendelianFilterInput,
    PlinkValidateInput,
)
from publicgenomicagent.tools.call_variants import call_variants
from publicgenomicagent.tools.mendelian import (
    mendelian_filter,
    plink_validate,
)


FIXTURES = Path(__file__).parents[1] / "fixtures"
MAKE_COHORT = FIXTURES / "make_cohort_fixture.py"
MAKE_BAM = FIXTURES / "make_patient_bam.py"
SAMTOOLS = str(Path("~/.pga/envs/pga-hts/bin/samtools").expanduser())
BCFTOOLS = str(Path("~/.pga/envs/pga-hts/bin/bcftools").expanduser())
REF_CACHE = Path("~/.pga/cache/reference").expanduser()
CONTIG = "chr2_roi"


def _env_ready() -> bool:
    try:
        reg = load_registry(registry_path())
        runtime = ToolRuntime(reg)
        runtime.version("plink")
        runtime.version("bcftools")
        return True
    except Exception:
        return False


def _count_variants(vcf: Path) -> int:
    if not vcf.exists():
        return 0
    n = 0
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                n += 1
    return n


def _query_variant(vcf: Path, pos: int) -> list[str]:
    """Devuelve las líneas del VCF en esa posición."""
    matches = []
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) > 1 and fields[1] == str(pos):
                matches.append(line.rstrip("\n"))
    return matches


@pytest.fixture
def trio_dir(tmp_path: Path) -> dict:
    """Genera el trío (C2, F1, M1) con cohort, referencia y VCFs."""
    work = tmp_path / "cohort"
    work.mkdir()

    # 1. Fixture
    subprocess.run(
        [sys.executable, str(MAKE_COHORT), str(work)],
        check=True, capture_output=True,
    )

    # 2. Referencia
    src = REF_CACHE / "hg38_chr2_110000001_110025000.fa"
    if not src.exists():
        pytest.skip(f"referencia no cacheada: {src}")
    ref = work / "nphp1_ref.fa"
    lines = src.read_text().splitlines()
    seq = "".join(l for l in lines if not l.startswith(">"))
    with ref.open("w") as f:
        f.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
    subprocess.run([SAMTOOLS, "faidx", str(ref)], check=True)

    # 3. BAMs del trío
    for s in ("C2", "F1", "M1"):
        sl = s.lower()
        subprocess.run(
            [sys.executable, str(MAKE_BAM),
             str(work / "cohort.vcf.gz"), str(ref), s,
             str(work / f"{sl}.bam"), "30"],
            check=True, capture_output=True,
        )
        # FASTQ no hace falta para el calling lineal

    # 4. Llamar variantes por individuo (lineal) usando la función
    # Python directamente (evita depender del CLI en PATH).
    reg = load_registry(registry_path())
    runtime = ToolRuntime(reg)
    for s in ("C2", "F1", "M1"):
        sl = s.lower()
        call_variants(runtime, CallVariantsInput(
            bam_path=work / f"{sl}.bam",
            reference_fasta=ref,
            output_vcf=work / f"{sl}.vcf.gz",
            region=f"{CONTIG}:1-25000",
            min_qual=0, min_dp=1,
        ))

    # 5. Renombrar sample y merge
    for s in ("C2", "F1", "M1"):
        sl = s.lower()
        # Extraer y renombrar
        result = subprocess.run(
            [BCFTOOLS, "view", str(work / f"{sl}.vcf.gz")],
            check=True, capture_output=True, text=True,
        )
        text = result.stdout.replace("\tSAMPLE", f"\t{s}")
        named = work / f"{sl}_named.vcf"
        named.write_text(text)
        subprocess.run(
            [BCFTOOLS, "view", str(named), "-Oz", "-o", str(work / f"{sl}_named.vcf.gz")],
            check=True, capture_output=True,
        )
        subprocess.run(
            [BCFTOOLS, "index", "-t", str(work / f"{sl}_named.vcf.gz")],
            check=True, capture_output=True,
        )

    subprocess.run(
        [BCFTOOLS, "merge", "-m", "all",
         str(work / "c2_named.vcf.gz"),
         str(work / "f1_named.vcf.gz"),
         str(work / "m1_named.vcf.gz"),
         "-Oz", "-o", str(work / "trio.vcf.gz")],
        check=True, capture_output=True,
    )
    subprocess.run(
        [BCFTOOLS, "index", "-t", str(work / "trio.vcf.gz")],
        check=True, capture_output=True,
    )

    # 6. .fam
    fam = work / "trio.fam"
    fam.write_text(
        "F1\tF1\t0\t0\t1\t1\n"
        "M1\tM1\t0\t0\t2\t1\n"
        "C2\tC2\tF1\tM1\t1\t2\n"
    )

    return {
        "tmp": tmp_path,
        "trio_vcf": work / "trio.vcf.gz",
        "fam": fam,
    }


def _pedigree() -> list[IndividualSpec]:
    return [
        IndividualSpec(sample="F1", sex="M", affected=False),
        IndividualSpec(sample="M1", sex="F", affected=False),
        IndividualSpec(sample="C2", sex="M", affected=True,
                       father="F1", mother="M1"),
    ]


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-mendelian y pga-hts")
def test_mendelian_filter_de_novo_zero(trio_dir: dict) -> None:
    """No hay variantes de novo en el fixture (todos los padres
    portan o no portan, sin hijo único)."""
    out = trio_dir["tmp"] / "mendel_out"
    result = mendelian_filter(MendelianFilterInput(
        trio_vcf=trio_dir["trio_vcf"],
        output_dir=out,
        pedigree=_pedigree(),
        proband="C2",
        min_dp=5, min_qual=0,
    ))
    assert result.counts["de_novo"] == 0


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-mendelian y pga-hts")
def test_mendelian_filter_auto_dom_relaxed_includes_target(trio_dir: dict) -> None:
    """Con filter_by_affected=False, la variante 8979 (padre portador
    sano, hijo portador afectado) aparece en auto_dom."""
    out = trio_dir["tmp"] / "mendel_relaxed"
    result = mendelian_filter(MendelianFilterInput(
        trio_vcf=trio_dir["trio_vcf"],
        output_dir=out,
        pedigree=_pedigree(),
        proband="C2",
        min_dp=5, min_qual=0,
        filter_by_affected=False,
    ))
    assert result.counts["auto_dom"] > 0
    # 8979 debe estar
    matches = _query_variant(result.auto_dom_vcf, 8979)
    assert len(matches) == 1
    fields = matches[0].split("\t")
    assert fields[3] == "A" and fields[4] == "T"


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-mendelian y pga-hts")
def test_mendelian_filter_auto_dom_strict_excludes_target(trio_dir: dict) -> None:
    """Con filter_by_affected=True, la variante 8979 NO aparece
    porque el padre portador está sano."""
    out = trio_dir["tmp"] / "mendel_strict"
    result = mendelian_filter(MendelianFilterInput(
        trio_vcf=trio_dir["trio_vcf"],
        output_dir=out,
        pedigree=_pedigree(),
        proband="C2",
        min_dp=5, min_qual=0,
        filter_by_affected=True,
    ))
    # 8979 no debe estar
    matches = _query_variant(result.auto_dom_vcf, 8979)
    assert len(matches) == 0


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-mendelian y pga-hts")
def test_mendelian_filter_high_dp_filters_all(trio_dir: dict) -> None:
    """Con min_dp=1000, ninguna variante pasa el filtro."""
    out = trio_dir["tmp"] / "mendel_high_dp"
    result = mendelian_filter(MendelianFilterInput(
        trio_vcf=trio_dir["trio_vcf"],
        output_dir=out,
        pedigree=_pedigree(),
        proband="C2",
        min_dp=1000, min_qual=0,
        filter_by_affected=False,
    ))
    assert result.counts["auto_dom"] == 0
    assert result.counts["de_novo"] == 0
    assert result.counts["auto_rec_hom"] == 0


@pytest.mark.integration
@pytest.mark.skipif(not _env_ready(), reason="requiere pga-mendelian y pga-hts")
def test_plink_validate_runs(trio_dir: dict) -> None:
    """PLINK convierte el VCF, corre --mendel sin errores y --genome."""
    out = trio_dir["tmp"] / "plink_out"
    result = plink_validate(PlinkValidateInput(
        trio_vcf=trio_dir["trio_vcf"],
        output_dir=out,
        pedigree=_pedigree(),
        proband="C2",
    ))
    assert result.ok
    assert result.bed_prefix.exists() or Path(str(result.bed_prefix) + ".bed").exists()
    assert result.n_mendel_errors == 0
    assert result.ibd_report.exists()
