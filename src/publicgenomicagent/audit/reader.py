"""Carga de un session.json en un AgentState y metadatos asociados."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from ..agent.state import AgentState


@dataclass(slots=True)
class Session:
    """Sesión cargada: estado + ruta de origen."""

    state: AgentState
    path: Path


def load_session(path: Path) -> Session:
    """Carga un session.json y lo valida contra AgentState.

    Levanta FileNotFoundError si no existe, ValueError si el JSON
    no valida contra el modelo, o json.JSONDecodeError si no es JSON.
    """
    if not path.exists():
        raise FileNotFoundError(f"no existe: {path}")

    raw = json.loads(path.read_text(encoding="utf-8"))
    state = AgentState.model_validate(raw)
    return Session(state=state, path=path)
