"""Puerta de calidad: ningún análisis sin verificación de entrada.

Principio: primero QC y verificación, después análisis y representación.
La puerta está ABIERTA solo si cada BAM del caso tiene una llamada
`qc_bam` correcta (ok) cuyo output dice `passed: True`. Cualquier otra
situación la deja CERRADA:

  - pendiente: un BAM sin QC (o cuyo QC se cayó);
  - fallida:   un BAM con `passed: False` (o sin el campo `passed`).

Un caso sin BAMs (solo fenotipo) no tiene nada que verificar: abierta.

La comprueban dos sitios:
  - los planners, para no proponer análisis con la puerta cerrada;
  - `AgentLoop`, como único punto por el que pasa toda ejecución
    (planner, LLM o Steps clásicos): una tool que no sea de QC no se
    ejecuta si la puerta no está abierta.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..tools.base import QCLevel
from .state import AgentState

# Tools que pertenecen a la capa de verificación (pueden correr siempre).
QC_TOOLS: frozenset[str] = frozenset({"qc_bam"})


@dataclass
class GateStatus:
    open: bool
    pending: list[str] = field(default_factory=list)          # samples sin QC
    failed: dict[str, list[str]] = field(default_factory=dict)  # sample -> issues

    def reason(self) -> str:
        parts = []
        if self.failed:
            parts.append("QC fallido: " + "; ".join(
                f"{s}: {', '.join(i) or 'sin detalle'}" for s, i in self.failed.items()))
        if self.pending:
            parts.append(f"QC pendiente: {self.pending}")
        return " | ".join(parts) or "QC aprobado"


def _samples_with_bams(state: AgentState) -> list[tuple[str, Path]]:
    if state.case.pedigree:
        names = [ind.sample for ind in state.case.pedigree]
    else:
        names = list(state.bams.keys())
    return [(s, state.bams[s]) for s in names if s in state.bams]


def qc_gate(state: AgentState) -> GateStatus:
    """Estado de la puerta para el estado actual de la sesión."""
    pending: list[str] = []
    failed: dict[str, list[str]] = {}
    for sample, bam in _samples_with_bams(state):
        calls = [c for c in state.calls_of("qc_bam")
                 if c.ok and c.input.get("bam_path") == str(bam)]
        if not calls:
            pending.append(sample)
            continue
        out = calls[-1].output or {}
        if out.get("passed") is not True:
            failed[sample] = list(out.get("issues") or ["output sin passed=True"])
    return GateStatus(open=not pending and not failed, pending=pending, failed=failed)


def qc_args(state: AgentState, sample: str, bam: Path) -> dict[str, Any]:
    """Argumentos de `qc_bam` para la puerta: recuentos + muestra + referencia."""
    args: dict[str, Any] = {
        "bam_path": bam,
        "level": QCLevel.COUNTS.value,
        "expected_sample": sample,
    }
    ref = state.references.get("hg38") or state.references.get("ref")
    if ref is not None:
        fai = Path(str(ref) + ".fai")
        if fai.exists():
            args["expected_reference_fai"] = fai
    return args
