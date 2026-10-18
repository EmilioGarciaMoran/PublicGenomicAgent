"""Validaciones de integridad de una sesión.

Devuelve una lista de problemas (vacía si todo OK). No levanta
excepciones: el caller decide qué hacer con los issues.

Incluye:
  - Verificación de hashes SHA-256 de inputs/outputs (audit v2).
  - Monotonía de timestamps.
  - Coherencia ok/output/error.
  - Validación de inputs contra el modelo Pydantic de cada tool.
"""
from __future__ import annotations

import hashlib
import json

from ..tools.registry import TOOL_REGISTRY, get_tool
from .reader import Session


def _json_canonical(obj) -> str:
    """Serialización canónica para hashing (idéntica a state.py)."""
    return json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def _sha256(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def verify_session(session: Session) -> list[str]:
    """Devuelve la lista de problemas detectados en la sesión."""
    issues: list[str] = []
    st = session.state

    # 0. Verificación de integridad por hashes (audit v2).
    # Solo se verifica si el record tiene hashes (sesiones antiguas no).
    for i, c in enumerate(st.tool_calls):
        if c.input_hash is not None:
            recomputed = _sha256(_json_canonical(c.input))
            if recomputed != c.input_hash:
                issues.append(
                    f"tool_call #{i} ({c.tool_name}): input_hash no coincide. "
                    f"Esperado {c.input_hash[:16]}..., "
                    f"recalculado {recomputed[:16]}... "
                    f"(posible manipulación del session.json)"
                )
        if c.output_hash is not None and c.output is not None:
            recomputed = _sha256(_json_canonical(c.output))
            if recomputed != c.output_hash:
                issues.append(
                    f"tool_call #{i} ({c.tool_name}): output_hash no coincide. "
                    f"Esperado {c.output_hash[:16]}..., "
                    f"recalculado {recomputed[:16]}... "
                    f"(posible manipulación del session.json)"
                )

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

