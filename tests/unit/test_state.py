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


# ---------------------------------------------------------------------------
# Privacidad: las claves de vcfs/artifacts no contienen rutas absolutas
# ---------------------------------------------------------------------------

def test_side_effects_use_logical_labels(tmp_path: Path):
    """Los efectos colaterales indexan etiquetas, no rutas."""
    state = AgentState(case=_minimal_case())
    ctx = SessionContext(state=state, runtime=None)

    # Simulamos el efecto colateral de un fetch_roi y un call_variants
    from publicgenomicagent.tools.base import FetchROIOutput, CallVariantsOutput
    from publicgenomicagent.agent.state import _register_side_effects

    roi_out = FetchROIOutput(
        ok=True,
        tool="fetch_roi",
        message="",
        outputs={},
        output_bam=tmp_path / "father.chr1_1000_1200.bam",
        output_bai=tmp_path / "father.chr1_1000_1200.bam.bai",
        bytes_written=123,
    )
    _register_side_effects(state, "fetch_roi", roi_out)

    cv_out = CallVariantsOutput(
        ok=True,
        tool="call_variants",
        message="",
        outputs={},
        output_vcf=tmp_path / "father.chr1_1000_1200.vcf.gz",
        output_tbi=tmp_path / "father.chr1_1000_1200.vcf.gz.tbi",
        variants_total=1,
        variants_passing=1,
    )
    _register_side_effects(state, "call_variants", cv_out)

    # Las claves deben ser etiquetas lógicas, no rutas
    artifact_keys = list(state.artifacts.keys())
    vcf_keys = list(state.vcfs.keys())

    assert any("roi_bam:father" == k for k in artifact_keys), artifact_keys
    assert any("called:father" == k for k in vcf_keys), vcf_keys

    # Y ninguna clave debe contener "/"
    for k in artifact_keys + vcf_keys:
        assert "/" not in k, f"clave con ruta absoluta: {k}"


def test_render_state_sanitizes_keys(tmp_path: Path):
    """Si una clave contiene una ruta (por un bug futuro), render_state la limpia."""
    from publicgenomicagent.agent.prompts import render_state

    state = AgentState(case=_minimal_case())
    # Inyectamos una clave con ruta absoluta directamente
    state.artifacts["roi_bam:/home/hospital/paciente_1234.exoma.bam"] = (
        tmp_path / "x.bam"
    )

    out = render_state(state, include_paths=False)

    # El prompt NO debe contener la ruta del hospital
    assert "/home/hospital/" not in out
    # Pero sí el basename
    assert "paciente_1234.exoma.bam" in out
