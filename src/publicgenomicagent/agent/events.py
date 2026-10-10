"""Eventos emitidos por el loop agéntico durante la ejecución.

Independientes de cualquier UI. `ui/` los traduce a render.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal
import time

EventKind = Literal[
    "step_start",
    "tool_call",
    "tool_result",
    "note",
    "final",
    "error",
]


@dataclass(slots=True)
class Event:
    kind: EventKind
    payload: dict[str, Any] = field(default_factory=dict)
    ts: float = field(default_factory=time.time)
