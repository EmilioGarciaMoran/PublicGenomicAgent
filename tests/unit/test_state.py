"""Tests del estado del agente y del SessionContext.

No requieren BAMs ni runtime real: usan tools que no necesitan
runtime (compare_vcfs, extract_hpo) para probar el registro de
trazas, y verifica que las tools que sí lo requieren fallan
de forma controlada cuando no hay runtime.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from publicgenomicagent.agent.state import (
    AgentState,
    SessionContext,
    ToolCallRecord,
)
from publicgenomicagent.tools.base import CaseManifest, GenomicRange


def _minimal_case() -> CaseManifest:
    return CaseManifest(
        case_id="TEST-001",
        source="clinical",
        candidate_rois=[GenomicRange(chrom="chr1", start=1, end=100, label="TEST")],
    )


def test_agent_state_minimal():
    state = AgentState(case=_minimal_case())
    assert state.case.case_id == "TEST-001"
    assert state.tool_calls == []
    assert state.bams == {}
    assert state.latest("fetch_roi") is None


def test_record_appends_and_indexes():
    state = AgentState(case=_minimal_case())
    rec = state.record(
        tool_name="fetch_roi",
        inp={"bam_path": "/a.bam", "region": "chr1:1-100", "output_bam": "/o.bam"},
        out={"output_bam": "/o.bam", "bytes_written": 42},
        duration_ms=15,
    )
    assert isinstance(rec, ToolCallRecord)
    assert rec.ok is True
    assert rec.error is None
    assert len(state.tool_calls) == 1
    assert state.latest("fetch_roi") == {
        "output_bam": "/o.bam",
        "bytes_written": 42,
    }


def test_record_error_does_not_index_output():
    state = AgentState(case=_minimal_case())
    rec = state.record(
        tool_name="fetch_roi",
        inp={"bam_path": "/a.bam"},
        out=None,
        duration_ms=5,
        error="RuntimeError: boom",
    )
    assert rec.ok is False
    assert rec.error == "RuntimeError: boom"
    assert state.latest("fetch_roi") is None


def test_calls_of_filters_by_tool_name():
    state = AgentState(case=_minimal_case())
    state.record(tool_name="a", inp={}, out={}, duration_ms=1)
    state.record(tool_name="b", inp={}, out={}, duration_ms=1)
    state.record(tool_name="a", inp={}, out={}, duration_ms=1)
    assert len(state.calls_of("a")) == 2
    assert len(state.calls_of("b")) == 1
    assert state.calls_of("z") == []


def test_has_tool():
    state = AgentState(case=_minimal_case())
    assert state.has_tool("x") is False
    state.record(tool_name="x", inp={}, out={"k": 1}, duration_ms=1)
    assert state.has_tool("x") is True


def test_note_appends():
    state = AgentState(case=_minimal_case())
    state.note("hola")
    state.note("mundo")
    assert state.notes == ["hola", "mundo"]


def test_to_json_roundtrip(tmp_path: Path):
    state = AgentState(case=_minimal_case())
    state.record(
        tool_name="fetch_roi",
        inp={"region": "chr1:1-100"},
        out={"bytes_written": 1},
        duration_ms=3,
    )
    out_path = tmp_path / "session.json"
    state.to_json(out_path)
    assert out_path.exists()
    data = json.loads(out_path.read_text())
    assert data["case"]["case_id"] == "TEST-001"
    assert len(data["tool_calls"]) == 1
    assert data["tool_calls"][0]["tool_name"] == "fetch_roi"


def test_session_context_call_tool_without_runtime_needed():
    """compare_vcfs no necesita runtime; debe ejecutarse sin runtime.

    Pero como compare_vcfs intenta abrir VCFs reales, forzamos
    un input inválido para que la tool falle y se registre el error.
    Así probamos el registro sin depender de ficheros.
    """
    state = AgentState(case=_minimal_case())
    ctx = SessionContext(state=state, runtime=None)

    with pytest.raises(Exception):
        ctx.call_tool(
            "compare_vcfs",
            baseline_vcf="/no/existe_a.vcf.gz",
            candidate_vcf="/no/existe_b.vcf.gz",
            output_delta_vcf="/tmp/delta.vcf.gz",
        )

    assert len(state.tool_calls) == 1
    rec = state.tool_calls[0]
    assert rec.tool_name == "compare_vcfs"
    assert rec.ok is False
    assert rec.error is not None
    assert "baseline_vcf" in rec.input


def test_session_context_requires_runtime_for_runtime_tools():
    state = AgentState(case=_minimal_case())
    ctx = SessionContext(state=state, runtime=None)

    with pytest.raises(ValueError, match="requiere un ToolRuntime"):
        ctx.call_tool(
            "fetch_roi",
            bam_path="/tmp/x.bam",
            region="chr1:1-100",
            output_bam="/tmp/x.roi.bam",
        )

    assert len(state.tool_calls) == 1
    rec = state.tool_calls[0]
    assert rec.tool_name == "fetch_roi"
    assert rec.ok is False
    assert "ToolRuntime" in (rec.error or "")


def test_session_context_invalid_input_is_recorded():
    state = AgentState(case=_minimal_case())
    ctx = SessionContext(state=state, runtime=None)

    with pytest.raises(Exception):
        ctx.call_tool("fetch_roi", region="chr1:1-100")  # faltan campos

    assert len(state.tool_calls) == 1
    rec = state.tool_calls[0]
    assert rec.ok is False
    assert "input validation failed" in (rec.error or "")
