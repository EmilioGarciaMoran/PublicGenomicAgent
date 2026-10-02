"""Tests de humo del registro de tools.

No requieren BAMs reales ni entornos micromamba. Solo validan:
  - que el registro se importa sin errores
  - que todas las tools declaradas existen
  - que dispatch valida inputs correctamente
  - que list_tools / describe_tool funcionan
"""
from __future__ import annotations

import pytest
from pydantic import ValidationError

from publicgenomicagent.tools.registry import (
    TOOL_REGISTRY,
    describe_tool,
    dispatch,
    get_tool,
    list_tools,
)

EXPECTED_TOOLS = {
    "fetch_roi",
    "qc_bam",
    "call_variants",
    "joint_call",
    "annotate_variants",
    "compare_vcfs",
    "local_pangenome",
    "build_local_graph",
    "align_to_graph",
    "call_from_graph",
    "mendelian_filter",
    "plink_validate",
    "extract_hpo",
    "phenotype_ranking",
}


def test_registry_imports():
    assert len(TOOL_REGISTRY) == len(EXPECTED_TOOLS)


def test_all_expected_tools_present():
    assert set(list_tools()) == EXPECTED_TOOLS


def test_list_tools_is_sorted():
    names = list_tools()
    assert names == sorted(names)


def test_get_tool_unknown_raises():
    with pytest.raises(KeyError, match="no registrada"):
        get_tool("no_existe")


def test_describe_tool_returns_schemas():
    d = describe_tool("fetch_roi")
    assert d["name"] == "fetch_roi"
    assert d["needs_runtime"] is True
    assert "properties" in d["input_schema"]
    assert "region" in d["input_schema"]["properties"]


def test_dispatch_unknown_tool_raises():
    with pytest.raises(KeyError):
        dispatch("no_existe", region="chr1:1-100")


def test_dispatch_invalid_input_raises_validation_error():
    # fetch_roi requiere bam_path, region, output_bam
    with pytest.raises(ValidationError):
        dispatch("fetch_roi", region="chr1:1-100")  # faltan bam_path y output_bam


def test_dispatch_validate_only_returns_input():
    inp = dispatch(
        "fetch_roi",
        validate_only=True,
        bam_path="/tmp/x.bam",
        region="chr1:1-100",
        output_bam="/tmp/x.roi.bam",
    )
    assert inp.region == "chr1:1-100"
    assert str(inp.bam_path) == "/tmp/x.bam"


def test_dispatch_tool_without_runtime_does_not_require_it():
    # compare_vcfs no necesita runtime; validate_only evita ejecutarla
    inp = dispatch(
        "compare_vcfs",
        validate_only=True,
        baseline_vcf="/tmp/a.vcf.gz",
        candidate_vcf="/tmp/b.vcf.gz",
        output_delta_vcf="/tmp/delta.vcf.gz",
    )
    assert str(inp.baseline_vcf) == "/tmp/a.vcf.gz"


def test_dispatch_requires_runtime_when_tool_needs_it():
    # fetch_roi necesita runtime; sin validate_only y sin runtime -> ValueError
    with pytest.raises(ValueError, match="requiere un ToolRuntime"):
        dispatch(
            "fetch_roi",
            bam_path="/tmp/x.bam",
            region="chr1:1-100",
            output_bam="/tmp/x.roi.bam",
        )


def test_all_tools_declare_nonempty_description():
    for name, spec in TOOL_REGISTRY.items():
        assert spec.description, f"{name} sin descripción"
        assert spec.tags, f"{name} sin tags"


def test_all_tools_have_pydantic_models():
    from pydantic import BaseModel
    for name, spec in TOOL_REGISTRY.items():
        assert issubclass(spec.input_model, BaseModel), name
        assert issubclass(spec.output_model, BaseModel), name
