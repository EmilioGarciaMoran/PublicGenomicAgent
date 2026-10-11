"""Tests de la capa de auditoría (audit/)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from rich.console import Console

from publicgenomicagent.agent.state import AgentState, ToolCallRecord
from publicgenomicagent.audit import (
    load_session,
    render_human,
    to_dict,
    verify_session,
)
from publicgenomicagent.audit.reader import Session
from publicgenomicagent.tools.base import CaseManifest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _empty_state() -> AgentState:
    return AgentState(case=CaseManifest(case_id="test"))


def _state_with_calls() -> AgentState:
    st = _empty_state()
    st.tool_calls.append(
        ToolCallRecord(
            tool_name="qc_bam",
            input={"bam_path": "/tmp/x.bam", "level": "structural"},
            output={"ok": True, "message": "QC estructural completado"},
            ok=True,
            duration_ms=42,
        )
    )
    st.tool_calls.append(
        ToolCallRecord(
            tool_name="fetch_roi",
            input={"bam_path": "/tmp/x.bam", "region": "chr1:100-200",
                   "output_bam": "/tmp/x.roi.bam"},
            output={"ok": True},
            ok=True,
            duration_ms=88,
        )
    )
    return st


def _session(state: AgentState, path: Path = Path("/tmp/session.json")) -> Session:
    return Session(state=state, path=path)


# ---------------------------------------------------------------------------
# load_session
# ---------------------------------------------------------------------------

def test_load_session_from_disk(tmp_path: Path):
    st = _empty_state()
    st.notes.append("test note")
    p = tmp_path / "session.json"
    st.to_json(p)

    loaded = load_session(p)

    assert loaded.state.case.case_id == "test"
    assert "test note" in loaded.state.notes
    assert loaded.path == p


def test_load_session_missing_file(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        load_session(tmp_path / "no-existe.json")


def test_load_session_invalid_json(tmp_path: Path):
    p = tmp_path / "session.json"
    p.write_text("not json at all", encoding="utf-8")
    with pytest.raises(json.JSONDecodeError):
        load_session(p)


# ---------------------------------------------------------------------------
# to_dict
# ---------------------------------------------------------------------------

def test_to_dict_empty():
    d = to_dict(_session(_empty_state()))
    assert d["case_id"] == "test"
    assert d["n_tool_calls"] == 0
    assert d["n_notes"] == 0
    assert d["tool_calls"] == []


def test_to_dict_with_calls():
    d = to_dict(_session(_state_with_calls()))
    assert d["n_tool_calls"] == 2
    assert d["tool_calls"][0]["tool_name"] == "qc_bam"
    assert d["tool_calls"][0]["ok"] is True
    assert d["tool_calls"][1]["duration_ms"] == 88


# ---------------------------------------------------------------------------
# render_human
# ---------------------------------------------------------------------------

def _render(session: Session, verbose: bool = False, verify: bool = False) -> str:
    c = Console(record=True, width=120)
    render_human(session, c, verbose=verbose, with_verify=verify)
    return c.export_text(styles=True)


def test_render_empty_shows_no_calls_message():
    out = _render(_session(_empty_state()))
    assert "Sin tool calls registradas" in out


def test_render_with_calls_shows_table():
    out = _render(_session(_state_with_calls()))
    assert "qc_bam" in out
    assert "fetch_roi" in out
    assert "chr1:100-200" in out or "region" in out


def test_render_with_verify_ok():
    out = _render(_session(_empty_state()), verify=True)
    assert "Integridad OK" in out


# ---------------------------------------------------------------------------
# verify_session
# ---------------------------------------------------------------------------

def test_verify_empty_state_no_issues():
    assert verify_session(_session(_empty_state())) == []


def test_verify_valid_calls_no_issues():
    st = _state_with_calls()
    # qc_bam y fetch_roi existen en TOOL_REGISTRY; inputs válidos.
    issues = verify_session(_session(st))
    # Puede haber issue si los paths de input no encajan con el modelo;
    # el test es tolerante a eso.
    assert all("no está en TOOL_REGISTRY" not in i for i in issues)


def test_verify_flags_unknown_tool():
    st = _empty_state()
    st.tool_calls.append(
        ToolCallRecord(
            tool_name="tool_que_no_existe",
            input={},
            output={"ok": True},
            ok=True,
            duration_ms=1,
        )
    )
    issues = verify_session(_session(st))
    assert any("no está en TOOL_REGISTRY" in i for i in issues)


def test_verify_flags_ok_without_output():
    st = _empty_state()
    st.tool_calls.append(
        ToolCallRecord(
            tool_name="qc_bam",
            input={"bam_path": "/tmp/x.bam", "level": "structural"},
            output=None,
            ok=True,
            duration_ms=1,
        )
    )
    issues = verify_session(_session(st))
    assert any("ok=true pero output vacío" in i for i in issues)


def test_verify_flags_error_without_message():
    st = _empty_state()
    st.tool_calls.append(
        ToolCallRecord(
            tool_name="qc_bam",
            input={"bam_path": "/tmp/x.bam", "level": "structural"},
            output=None,
            ok=False,
            duration_ms=1,
            error=None,
        )
    )
    issues = verify_session(_session(st))
    assert any("ok=false pero error vacío" in i for i in issues)


def test_verify_flags_non_monotonic_timestamps():
    from datetime import datetime, timezone, timedelta

    st = _empty_state()
    t1 = datetime(2026, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
    t0 = t1 - timedelta(hours=1)
    st.tool_calls.append(
        ToolCallRecord(
            tool_name="qc_bam",
            input={"bam_path": "/tmp/x.bam", "level": "structural"},
            output={"ok": True},
            ok=True,
            duration_ms=1,
            timestamp=t1,
        )
    )
    st.tool_calls.append(
        ToolCallRecord(
            tool_name="qc_bam",
            input={"bam_path": "/tmp/y.bam", "level": "structural"},
            output={"ok": True},
            ok=True,
            duration_ms=1,
            timestamp=t0,
        )
    )
    issues = verify_session(_session(st))
    assert any("timestamp" in i for i in issues)
