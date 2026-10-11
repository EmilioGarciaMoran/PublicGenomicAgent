"""Validaciones de integridad de una sesión.

Devuelve una lista de problemas (vacía si todo OK). No levanta
excepciones: el caller decide qué hacer con los issues.
"""
from __future__ import annotations

from ..tools.registry import TOOL_REGISTRY, get_tool
from .reader import Session


def verify_session(session: Session) -> list[str]:
    """Devuelve la lista de problemas detectados en la sesión."""
    issues: list[str] = []
    st = session.state

    # 1. Monotonía de timestamps en tool_calls
    prev_ts = None
    for i, c in enumerate(st.tool_calls):
        if prev_ts is not None and c.timestamp < prev_ts:
            issues.append(
                f"tool_call #{i} ({c.tool_name}): timestamp {c.timestamp} "
                f"anterior al previo {prev_ts}"
            )
        prev_ts = c.timestamp

    # 2. Cada tool_call: tool conocida, output/error coherentes
    for i, c in enumerate(st.tool_calls):
        if c.tool_name not in TOOL_REGISTRY:
            issues.append(
                f"tool_call #{i}: tool '{c.tool_name}' no está en TOOL_REGISTRY"
            )
            continue

        if c.ok and not c.output:
            issues.append(
                f"tool_call #{i} ({c.tool_name}): ok=true pero output vacío"
            )
        if not c.ok and not c.error:
            issues.append(
                f"tool_call #{i} ({c.tool_name}): ok=false pero error vacío"
            )

        # 3. Validación de input contra el modelo Pydantic de la tool
        try:
            spec = get_tool(c.tool_name)
            spec.input_model(**c.input)
        except Exception as e:  # noqa: BLE001
            issues.append(
                f"tool_call #{i} ({c.tool_name}): input inválido "
                f"({type(e).__name__}: {e})"
            )

    return issues
