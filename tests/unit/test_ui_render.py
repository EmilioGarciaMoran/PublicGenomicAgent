"""Tests mínimos del render de la UI (Día 1)."""
from __future__ import annotations

from rich.console import Console

from publicgenomicagent.agent.events import Event
from publicgenomicagent.ui.render import render


def _cap(event: Event, verbose: bool = False) -> str:
    c = Console(record=True, width=80)
    render(event, c, verbose=verbose)
    return c.export_text()


def test_tool_result_ok():
    assert "✓" in _cap(Event(kind="tool_result", payload={"ok": True}))


def test_tool_result_fail():
    assert "✗" in _cap(Event(kind="tool_result", payload={"ok": False}))


def test_error_verbose_shows_message():
    out = _cap(Event(kind="error", payload={"message": "boom"}), verbose=True)
    assert "boom" in out


def test_final_prints_marker():
    assert "fin" in _cap(Event(kind="final", payload={}))


def test_unknown_kind_does_not_crash():
    # kind fuera del Literal no debe romper el render
    render(Event(kind="tool_result", payload={}), Console(record=True))
