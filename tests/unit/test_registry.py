"""Tests del registro de tools del agente."""
from __future__ import annotations

import pytest

from publicgenomicagent.tools.registry import (
    TOOL_REGISTRY,
    describe_tool,
    dispatch,
    get_tool,
    list_tools,
)


EXPECTED_TOOLS = {
    "align_to_graph",
    "annotate_variants",
    "build_local_graph",
    "call_from_graph",
    "call_variants",
    "compare_vcfs",
    "detect_homozygous_regions",
    "extract_hpo",
    "fetch_roi",
    "joint_call",
    "local_pangenome",
    "mendelian_filter",
    "nl_to_case_manifest",
    "phenotype_ranking",
    "plink_validate",
    "prioritize_variants",
    "qc_bam",
}


def test_registry_imports():
    assert len(TOOL_REGISTRY) == len(EXPECTED_TOOLS)


def test_all_expected_tools_present():
    assert set(list_tools()) == EXPECTED_TOOLS


def test_list_tools_sorted():
    names = list_tools()
    assert names == sorted(names)


def test_get_tool_returns_spec():
    spec = get_tool("qc_bam")
    assert spec.name == "qc_bam"
    assert spec.needs_runtime is True


def test_get_tool_unknown_raises():
    with pytest.raises(KeyError):
        get_tool("no_existe_esta_tool")


def test_describe_tool_shape():
    d = describe_tool("prioritize_variants")
    assert d["name"] == "prioritize_variants"
    assert "input_schema" in d
    assert "output_schema" in d
    assert d["needs_runtime"] is False


def test_detect_homozygous_regions_registered():
    spec = get_tool("detect_homozygous_regions")
    assert spec.needs_runtime is True
    assert "homozygosity" in spec.tags
    assert "consanguinity" in spec.tags


def test_detect_homozygous_regions_input_schema():
    spec = get_tool("detect_homozygous_regions")
    props = spec.input_schema()["properties"]
    assert "vcf" in props
    assert "output_dir" in props
    assert "samples" in props
    assert "min_length_kb" in props
    assert "min_snps" in props
    assert "min_quality" in props


def test_detect_homozygous_regions_output_schema():
    spec = get_tool("detect_homozygous_regions")
    props = spec.output_schema()["properties"]
    assert "regions" in props
    assert "samples_analyzed" in props
    assert "total_regions" in props
    assert "total_mb" in props
    assert "f_roh_proxy" in props
    assert "caveat" in props


def test_dispatch_validate_only():
    inp = dispatch(
        "detect_homozygous_regions",
        validate_only=True,
        vcf="/tmp/x.vcf.gz",
        output_dir="/tmp/out",
    )
    assert inp.vcf is not None
    assert inp.min_length_kb == 100


def test_dispatch_runtime_required():
    with pytest.raises(ValueError, match="ToolRuntime"):
        dispatch(
            "detect_homozygous_regions",
            runtime=None,
            vcf="/tmp/x.vcf.gz",
            output_dir="/tmp/out",
        )
