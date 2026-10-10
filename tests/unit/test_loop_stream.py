"""Tests de AgentLoop.stream() — equivalencia con run()."""
from __future__ import annotations

from publicgenomicagent.agent.events import Event
from publicgenomicagent.agent.loop import AgentLoop
from publicgenomicagent.agent.state import AgentState, SessionContext
from publicgenomicagent.tools.base import CaseManifest


def _empty_ctx() -> SessionContext:
    return SessionContext(state=AgentState(case=CaseManifest(case_id="t")))


def test_stream_emits_final_on_empty_state():
    loop = AgentLoop()
    events = list(loop.stream(_empty_ctx()))
    assert events, "stream() no emitió nada"
    assert all(isinstance(e, Event) for e in events)
    assert events[-1].kind == "final"


def test_stream_yields_event_objects():
    loop = AgentLoop()
    for e in loop.stream(_empty_ctx()):
        assert isinstance(e, Event)
        assert e.kind in {
            "step_start", "tool_call", "tool_result",
            "note", "final", "error",
        }
