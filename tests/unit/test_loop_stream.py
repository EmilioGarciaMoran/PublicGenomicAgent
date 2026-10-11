"""Tests de AgentLoop.stream() — equivalencia con run() y no duplicación."""
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


class _CountingPlanner:
    """Planner falso que cuenta llamadas y siempre devuelve None (stop)."""

    def __init__(self) -> None:
        self.calls = 0

    def next_action(self, state):  # noqa: ANN001
        self.calls += 1
        return None


def test_stream_calls_planner_once_per_iteration():
    """Regresión: stream() no debe llamar al planner dos veces por iteración.

    Antes del fix, stream() consultaba al planner y luego llamaba a
    _run_one, que volvía a consultarlo. Con planner que devuelve None,
    el resultado eran dos invocaciones por iteración (doble latencia,
    doble nota en state.notes).
    """
    planner = _CountingPlanner()
    loop = AgentLoop(planner=planner)
    ctx = _empty_ctx()

    events = list(loop.stream(ctx))

    assert planner.calls == 1, (
        f"planner llamado {planner.calls} veces, esperado 1"
    )
    assert events[-1].kind == "final"

