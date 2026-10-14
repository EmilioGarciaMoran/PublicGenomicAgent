"""Tests unitarios de `detect_homozygous_regions`.

No ejecutan `bcftools roh` de verdad: mockean el `ToolRuntime` para
devolver salidas predefinidas y verifican el parseo, los filtros y
la agregación.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import pytest

from publicgenomicagent.tools.base import (
    DetectHomozygousRegionsInput,
    DetectHomozygousRegionsOutput,
    HomozygousRegion,
)
from publicgenomicagent.tools.roh import (
    CAVEAT,
    _filter_regions,
    _parse_roh_output,
    detect_homozygous_regions,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Cabecera típica de bcftools roh -O r
_ROH_HEADER = (
    "# This file was produced by: bcftools roh(1.19+htslib-1.19.1)\n"
    "# The command line was:\tbcftools roh ...\n"
    "#\n"
    "# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t"
    "[6]Length (bp)\t[7]Number of markers\t[8]Quality\n"
)


def _roh_line(sample: str, chrom: str, start: int, end: int, n_markers: int, quality: float) -> str:
    length = end - start + 1
    return f"RG\t{sample}\t{chrom}\t{start}\t{end}\t{length}\t{n_markers}\t{quality}\n"


def _fake_runtime(outputs: dict[tuple[str, ...], str]) -> MagicMock:
    """Crea un ToolRuntime mockeado.

    `outputs` mapea tuplas de args (sin el primer 'roh') a stdout a devolver.
    """
    runtime = MagicMock()

    def fake_run(tool: str, args: list[str], **kwargs: Any) -> MagicMock:
        result = MagicMock()
        result.returncode = 0
        # Buscar el output por coincidencia de 'bcftools query -l' o
        # 'roh -s <sample>'.
        args_key = tuple(args)
        result.stdout = outputs.get(args_key, "")
        result.stderr = ""

        # Si es `query -l`, devolver samples
        # Si es `roh -s <sample>`, escribir el fichero -o
        if args and args[0] == "roh" and "-o" in args:
            out_idx = args.index("-o") + 1
            out_path = Path(args[out_idx])
            out_path.write_text(result.stdout or "")

        return result

    runtime.run.side_effect = fake_run
    return runtime


def _make_vcf(tmp_path: Path, name: str = "test.vcf.gz") -> Path:
    """Crea un fichero dummy como VCF. Solo importa que exista."""
    p = tmp_path / name
    p.write_text("dummy\n")
    return p


# ---------------------------------------------------------------------------
# _parse_roh_output
# ---------------------------------------------------------------------------

def test_parse_roh_output_empty() -> None:
    text = _ROH_HEADER
    regions = _parse_roh_output(text, "father")
    assert regions == []


def test_parse_roh_output_one_region() -> None:
    text = _ROH_HEADER + _roh_line("father", "chr2", 40000, 100000, 25, 30.5)
    regions = _parse_roh_output(text, "father")
    assert len(regions) == 1
    r = regions[0]
    assert r.sample == "father"
    assert r.chrom == "chr2"
    assert r.start == 40000
    assert r.end == 100000
    assert r.length_bp == 60001
    assert r.length_kb == 60.001
    assert r.n_markers == 25
    assert r.quality == 30.5


def test_parse_roh_output_multiple_regions() -> None:
    text = (
        _ROH_HEADER
        + _roh_line("father", "chr2", 40000, 100000, 25, 30.5)
        + _roh_line("father", "chr2", 200000, 300000, 15, 22.0)
    )
    regions = _parse_roh_output(text, "father")
    assert len(regions) == 2
    assert regions[1].start == 200000


def test_parse_roh_output_skips_malformed() -> None:
    text = _ROH_HEADER + "RG\tfather\tchr2\n"  # faltan campos
    regions = _parse_roh_output(text, "father")
    assert regions == []


# ---------------------------------------------------------------------------
# _filter_regions
# ---------------------------------------------------------------------------

def _make_region(**kwargs: Any) -> HomozygousRegion:
    defaults = dict(
        sample="x",
        chrom="chr2",
        start=0,
        end=50000,
        length_bp=50001,
        length_kb=50.001,
        n_markers=20,
        quality=30.0,
    )
    defaults.update(kwargs)
    return HomozygousRegion(**defaults)


def test_filter_by_min_length() -> None:
    short = _make_region(length_bp=10_000, length_kb=10.0)
    long = _make_region(length_bp=200_000, length_kb=200.0)
    kept = _filter_regions([short, long], min_length_kb=100, min_snps=1, min_quality=0.0)
    assert kept == [long]


def test_filter_by_min_snps() -> None:
    few = _make_region(n_markers=5)
    many = _make_region(n_markers=25)
    kept = _filter_regions([few, many], min_length_kb=0, min_snps=10, min_quality=0.0)
    assert kept == [many]


def test_filter_by_min_quality() -> None:
    low = _make_region(quality=5.0)
    high = _make_region(quality=20.0)
    kept = _filter_regions([low, high], min_length_kb=0, min_snps=0, min_quality=10.0)
    assert kept == [high]


def test_filter_all_out() -> None:
    r = _make_region(length_bp=5_000, length_kb=5.0, n_markers=3, quality=2.0)
    kept = _filter_regions([r], min_length_kb=100, min_snps=10, min_quality=10.0)
    assert kept == []


# ---------------------------------------------------------------------------
# detect_homozygous_regions (con runtime mockeado)
# ---------------------------------------------------------------------------

def test_detect_no_roh(tmp_path: Path) -> None:
    vcf = _make_vcf(tmp_path)
    out_dir = tmp_path / "out"

    runtime = _fake_runtime({
        ("query", "-l", str(vcf)): "father\n",
        ("roh", "-s", "father", "-e", "-", "-O", "r",
         "-o", str(out_dir / "roh.father.txt"), "-I", str(vcf)): _ROH_HEADER,
    })

    result = detect_homozygous_regions(
        DetectHomozygousRegionsInput(vcf=vcf, output_dir=out_dir),
        runtime,
    )

    assert isinstance(result, DetectHomozygousRegionsOutput)
    assert result.regions == []
    assert result.total_regions == 0
    assert result.total_mb == 0.0
    assert result.caveat == CAVEAT


def test_detect_with_one_region(tmp_path: Path) -> None:
    vcf = _make_vcf(tmp_path)
    out_dir = tmp_path / "out"

    stdout = _ROH_HEADER + _roh_line("father", "chr2", 40000, 200000, 30, 25.0)

    runtime = _fake_runtime({
        ("query", "-l", str(vcf)): "father\n",
        ("roh", "-s", "father", "-e", "-", "-O", "r",
         "-o", str(out_dir / "roh.father.txt"), "-I", str(vcf)): stdout,
    })

    result = detect_homozygous_regions(
        DetectHomozygousRegionsInput(vcf=vcf, output_dir=out_dir),
        runtime,
    )

    assert result.total_regions == 1
    assert result.regions[0].sample == "father"
    assert result.regions[0].length_bp == 160001
    assert result.total_mb > 0.0


def test_detect_multi_sample(tmp_path: Path) -> None:
    vcf = _make_vcf(tmp_path)
    out_dir = tmp_path / "out"

    f_out = _ROH_HEADER + _roh_line("father", "chr2", 40000, 200000, 30, 25.0)
    m_out = _ROH_HEADER  # sin ROH
    p_out = _ROH_HEADER + _roh_line("proband", "chr2", 50000, 250000, 40, 30.0)

    runtime = _fake_runtime({
        ("query", "-l", str(vcf)): "father\nmother\nproband\n",
        ("roh", "-s", "father", "-e", "-", "-O", "r",
         "-o", str(out_dir / "roh.father.txt"), "-I", str(vcf)): f_out,
        ("roh", "-s", "mother", "-e", "-", "-O", "r",
         "-o", str(out_dir / "roh.mother.txt"), "-I", str(vcf)): m_out,
        ("roh", "-s", "proband", "-e", "-", "-O", "r",
         "-o", str(out_dir / "roh.proband.txt"), "-I", str(vcf)): p_out,
    })

    result = detect_homozygous_regions(
        DetectHomozygousRegionsInput(vcf=vcf, output_dir=out_dir),
        runtime,
    )

    assert result.samples_analyzed == ["father", "mother", "proband"]
    assert result.total_regions == 2
    assert {r.sample for r in result.regions} == {"father", "proband"}


def test_detect_filters_short_regions(tmp_path: Path) -> None:
    vcf = _make_vcf(tmp_path)
    out_dir = tmp_path / "out"

    # Región de 50 kb (<100 kb por defecto)
    stdout = _ROH_HEADER + _roh_line("father", "chr2", 40000, 90000, 15, 20.0)

    runtime = _fake_runtime({
        ("query", "-l", str(vcf)): "father\n",
        ("roh", "-s", "father", "-e", "-", "-O", "r",
         "-o", str(out_dir / "roh.father.txt"), "-I", str(vcf)): stdout,
    })

    result = detect_homozygous_regions(
        DetectHomozygousRegionsInput(vcf=vcf, output_dir=out_dir),
        runtime,
    )

    assert result.total_regions == 0


def test_detect_vcf_not_found(tmp_path: Path) -> None:
    runtime = MagicMock()
    with pytest.raises(FileNotFoundError):
        detect_homozygous_regions(
            DetectHomozygousRegionsInput(
                vcf=tmp_path / "no-existe.vcf.gz",
                output_dir=tmp_path / "out",
            ),
            runtime,
        )


def test_detect_explicit_samples_skips_query(tmp_path: Path) -> None:
    vcf = _make_vcf(tmp_path)
    out_dir = tmp_path / "out"

    stdout = _ROH_HEADER
    runtime = _fake_runtime({
        ("roh", "-s", "father", "-e", "-", "-O", "r",
         "-o", str(out_dir / "roh.father.txt"), "-I", str(vcf)): stdout,
    })

    result = detect_homozygous_regions(
        DetectHomozygousRegionsInput(
            vcf=vcf, output_dir=out_dir, samples=["father"],
        ),
        runtime,
    )

    assert result.samples_analyzed == ["father"]
    # No se ha llamado a `bcftools query -l`
    for call in runtime.run.call_args_list:
        args = call[0][1] if len(call[0]) > 1 else call.kwargs.get("args", [])
        assert args[0] != "query"
