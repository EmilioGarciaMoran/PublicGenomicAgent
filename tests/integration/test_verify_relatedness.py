"""verify_relatedness con el bcftools real de pga-hts sobre BAMs simulados."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

from publicgenomicagent.env.micromamba import bin_path
from publicgenomicagent.env.paths import registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.tools.base import VerifyRelatednessInput
from publicgenomicagent.tools.pedigree import parse_ped_text
from publicgenomicagent.tools.relatedness import verify_relatedness

sys.path.insert(0, str(Path(__file__).parents[1] / "fixtures"))

pysam = pytest.importorskip("pysam")
pytestmark = pytest.mark.integration

TRIO_PED = "FAM F 0 0 1 1\nFAM M 0 0 2 1\nFAM P F M 1 2\n"
TRIO = {"F": (None, None), "M": (None, None), "P": ("F", "M"), "X": (None, None)}


def _runtime() -> ToolRuntime:
    if not bin_path("pga-hts", "bcftools").exists():
        pytest.skip("bcftools de pga-hts no disponible")
    return ToolRuntime(load_registry(registry_path()))


def _run(tmp_path: Path, bam_as=None, n_sites=300):
    from make_family_bams import make_family

    fam = make_family(tmp_path / "fam", TRIO, n_sites=n_sites, bam_as=bam_as)
    _, ped = parse_ped_text(TRIO_PED)
    inp = VerifyRelatednessInput(
        bams=[(s, fam["bams"][s]) for s in ("F", "M", "P")],
        reference_fasta=fam["ref"], sites=fam["sites"], pedigree=ped,
        output_dir=tmp_path / "out",
    )
    return verify_relatedness(_runtime(), inp)


def test_correct_trio_passes(tmp_path):
    out = _run(tmp_path)
    assert out.status == "passed", out.issues
    po = [p for p in out.pairs if p["declared"] == "padre/madre-hijo"]
    assert all(p["ibs0"] == 0 and p["n"] >= 250 for p in po)


def test_swapped_father_bam_fails_even_with_correct_sm_tag(tmp_path):
    out = _run(tmp_path, bam_as={"F": "X"})          # SM:F, lecturas de X
    assert out.status == "failed"
    assert any("incompatible con padre/madre-hijo" in i for i in out.issues)


def test_small_panel_is_not_verifiable(tmp_path):
    out = _run(tmp_path, n_sites=40)
    assert out.status == "not_verifiable"
