"""Tests for the Sandbox consumer module."""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from publicgenomicagent.agent.sandbox import (
    ZYGOSITY_TO_GT,
    SandboxCase,
    SandboxCatalog,
    SandboxError,
    _parse_spdi,
)


# --- helpers ----------------------------------------------------------

def _make_case(root: Path, gene: str, spdi: str, chrom: str,
               zygosity: dict[str, str], hpo: list[dict] | None = None,
               with_bams: bool = True) -> Path:
    """Create a minimal Sandbox case in `root/genes/<gene>/<spdi_norm>/`."""
    # Normalize SPDI to folder name (replace : with _)
    folder_name = spdi.replace(":", "_")
    case_dir = root / "genes" / gene / folder_name
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "gene": gene,
        "spdi": spdi,
        "title": f"Test case {gene}",
        "clinical_text": "Test clinical text.",
        "ground_truth": {
            "chrom": chrom,
            "zygosity": zygosity,
        },
        "hpo_expected": hpo or [],
    }
    (case_dir / "manifest.yaml").write_text(
        yaml.safe_dump(manifest), encoding="utf-8"
    )

    if with_bams:
        bams_dir = case_dir / "bams"
        bams_dir.mkdir(exist_ok=True)
        for sample in ("father", "mother", "proband"):
            (bams_dir / f"{sample}.bam").write_bytes(b"")
            (bams_dir / f"{sample}.bam.bai").write_bytes(b"")

    return case_dir


# --- SPDI parsing -----------------------------------------------------

def test_parse_spdi_snv():
    refseq, pos, del_len, ins = _parse_spdi("NC_000001.11:1050:T:G")
    assert refseq == "NC_000001.11"
    assert pos == 1050
    assert del_len == 1  # len("T")
    assert ins == "G"


def test_parse_spdi_deletion():
    refseq, pos, del_len, ins = _parse_spdi("NC_000002.12:110800000:290000:")
    assert refseq == "NC_000002.12"
    assert pos == 110800000
    assert del_len == 290000
    assert ins == ""


def test_parse_spdi_control():
    refseq, pos, del_len, ins = _parse_spdi("NC_000002.12:110800000::")
    assert refseq == "NC_000002.12"
    assert pos == 110800000
    assert del_len == 0
    assert ins == ""


def test_parse_spdi_invalid_raises():
    with pytest.raises(SandboxError):
        _parse_spdi("not-an-spdi")


# --- Zygosity mapping -------------------------------------------------

def test_zygosity_to_gt_mapping():
    assert ZYGOSITY_TO_GT["HOM_WT"] == "0/0"
    assert ZYGOSITY_TO_GT["HET"] == "0/1"
    assert ZYGOSITY_TO_GT["HOM_ALT"] == "1/1"
    assert ZYGOSITY_TO_GT["HOM_DEL"] == "1/1"


# --- SandboxCase methods ----------------------------------------------

def _dummy_case() -> SandboxCase:
    return SandboxCase(
        case_id="NPHP1:NC_000002.12:110800000:290000:",
        gene="NPHP1",
        spdi="NC_000002.12:110800000:290000:",
        title="Test",
        clinical_text="",
        chrom="chr2",
        zygosity={"father": "HET", "mother": "HET", "proband": "HOM_DEL"},
        hpo_expected=[],
        path=Path("/tmp"),
        bams={},
    )


def test_expected_genotypes():
    c = _dummy_case()
    assert c.expected_genotypes() == {"father": "0/1", "mother": "0/1", "proband": "1/1"}


def test_expected_proband_gt_hom_del():
    c = _dummy_case()
    assert c.expected_proband_gt() == "1/1"


def test_expected_proband_gt_hom_wt():
    c = SandboxCase(
        case_id="x", gene="x", spdi="x", title="", clinical_text="",
        chrom="chr1",
        zygosity={"father": "HOM_WT", "mother": "HOM_WT", "proband": "HOM_WT"},
        hpo_expected=[], path=Path("/tmp"), bams={},
    )
    assert c.expected_proband_gt() == "0/0"


def test_sandbox_zygosity_lookup():
    c = _dummy_case()
    assert c.sandbox_zygosity("father") == "HET"
    assert c.sandbox_zygosity("nope") is None


# --- SandboxCatalog ---------------------------------------------------

def test_catalog_from_filesystem(tmp_path: Path):
    _make_case(tmp_path, "NPHP1", "NC_000002.12:110800000:290000:",
               "chr2", {"father": "HET", "mother": "HET", "proband": "HOM_DEL"})
    _make_case(tmp_path, "CYP2D6", "NC_000022.11:42128940:C:T",
               "chr22", {"father": "HET", "mother": "HET", "proband": "HOM_ALT"})

    catalog = SandboxCatalog(tmp_path)
    cases = catalog.list_cases()
    assert len(cases) == 2
    genes = {c.gene for c in cases}
    assert genes == {"NPHP1", "CYP2D6"}


def test_catalog_get_by_gene(tmp_path: Path):
    _make_case(tmp_path, "NPHP1", "NC_000002.12:110800000:290000:",
               "chr2", {"father": "HET", "mother": "HET", "proband": "HOM_DEL"})
    catalog = SandboxCatalog(tmp_path)
    case = catalog.get_case("NPHP1")
    assert case.gene == "NPHP1"


def test_catalog_get_by_case_id(tmp_path: Path):
    _make_case(tmp_path, "NPHP1", "NC_000002.12:110800000:290000:",
               "chr2", {"father": "HET", "mother": "HET", "proband": "HOM_DEL"})
    catalog = SandboxCatalog(tmp_path)
    case = catalog.get_case("NPHP1:NC_000002.12:110800000:290000:")
    assert case.gene == "NPHP1"


def test_catalog_get_unknown_raises(tmp_path: Path):
    _make_case(tmp_path, "NPHP1", "NC_000002.12:110800000:290000:",
               "chr2", {"father": "HET", "mother": "HET", "proband": "HOM_DEL"})
    catalog = SandboxCatalog(tmp_path)
    with pytest.raises(SandboxError, match="No case found"):
        catalog.get_case("NOPE")


def test_catalog_missing_root_raises(tmp_path: Path):
    with pytest.raises(SandboxError, match="does not exist"):
        SandboxCatalog(tmp_path / "nope")


# --- to_case_manifest -------------------------------------------------

def test_to_case_manifest_trio(tmp_path: Path):
    _make_case(
        tmp_path, "NPHP1", "NC_000002.12:110800000:290000:", "chr2",
        {"father": "HET", "mother": "HET", "proband": "HOM_DEL"},
        hpo=[{"code": "HP:0003774", "term": "Chronic kidney disease"}],
    )
    catalog = SandboxCatalog(tmp_path)
    case = catalog.get_case("NPHP1")
    manifest = catalog.to_case_manifest(case)

    # Pedigree
    assert len(manifest.pedigree) == 3
    samples = [i.sample for i in manifest.pedigree]
    assert samples == ["father", "mother", "proband"]
    father = manifest.pedigree[0]
    proband = manifest.pedigree[2]
    assert father.affected is False
    assert proband.affected is True
    assert proband.father == "father"
    assert proband.mother == "mother"

    # ROI
    assert len(manifest.candidate_rois) == 1
    roi = manifest.candidate_rois[0]
    assert roi.chrom == "chr2"
    assert roi.start <= 110800000 <= roi.end
    assert roi.label == "NPHP1"

    # HPO
    assert "proband" in manifest.hpo_terms
    assert "HP:0003774" in manifest.hpo_terms["proband"]

    # case_id no debe contener ":" (válido como nombre de fichero)
    assert ":" not in manifest.case_id


def test_to_case_manifest_control_unaffected(tmp_path: Path):
    _make_case(
        tmp_path, "CONTROL", "NC_000002.12:110800000::", "chr2",
        {"father": "HOM_WT", "mother": "HOM_WT", "proband": "HOM_WT"},
    )
    catalog = SandboxCatalog(tmp_path)
    case = catalog.get_case("CONTROL")
    manifest = catalog.to_case_manifest(case)

    proband = manifest.pedigree[-1]
    assert proband.affected is False
    assert manifest.affected_samples == []

