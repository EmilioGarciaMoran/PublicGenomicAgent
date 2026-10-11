"""Auditoría de sesiones del agente.

Lee un `session.json` producido por `pga chat` (o por `pga run-trio`)
y lo convierte en un informe auditable: tabla de tools ejecutadas,
notas del loop, y validación de integridad opcional.

No modifica la sesión. Solo lectura.
"""
from .reader import load_session
from .report import render_human, to_dict
from .verify import verify_session

__all__ = ["load_session", "render_human", "to_dict", "verify_session"]
