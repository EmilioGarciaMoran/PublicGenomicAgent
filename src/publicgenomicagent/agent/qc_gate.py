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
QC_TOOLS: frozenset[str] = frozenset({"qc_bam", "verify_relatedness"})

# Excepción humana: aceptar una verificación familiar "no verificable".
OVERRIDE_UNVERIFIABLE = "relatedness_not_verifiable"


@dataclass
class GateStatus:
    open: bool
    pending: list[str] = field(default_factory=list)          # samples sin QC
    failed: dict[str, list[str]] = field(default_factory=dict)  # sample -> issues
    bams_ok: bool = False                                       # QC de todos los BAM aprobado
    relatedness: str | None = None    # None (no exigida) | pending | passed | warning | failed | not_verifiable

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


def relatedness_required(state: AgentState) -> bool:
    """Se exige verificación familiar si hay panel de sitios y 2+ secuenciados del PED."""
    if "relatedness_sites" not in state.references or not state.case.pedigree:
        return False
    return len(_samples_with_bams(state)) >= 2


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
    bams_ok = not pending and not failed

    rel: str | None = None
    if relatedness_required(state):
        calls = [c for c in state.calls_of("verify_relatedness") if c.ok]
        if not calls:
            rel = "pending"
            if bams_ok:
                pending.append("verificación familiar")
        else:
            out = calls[-1].output or {}
            rel = out.get("status") or "failed"
            if rel == "failed":
                failed["familia"] = list(out.get("issues") or ["verificación familiar fallida"])
            elif rel == "not_verifiable" and OVERRIDE_UNVERIFIABLE not in state.gate_overrides:
                failed["familia"] = ["verificación familiar no verificable (pocos sitios); "
                                     "requiere excepción humana registrada"]
            elif rel not in ("passed", "warning", "not_verifiable"):
                failed["familia"] = [f"estado desconocido: {rel}"]
    return GateStatus(open=not pending and not failed, pending=pending,
                      failed=failed, bams_ok=bams_ok, relatedness=rel)


def relatedness_args(state: AgentState) -> dict[str, Any]:
    """Argumentos de `verify_relatedness` para el caso actual."""
    ref = state.references.get("hg38") or state.references.get("ref")
    return {
        "bams": [(s, b) for s, b in _samples_with_bams(state)],
        "reference_fasta": ref,
        "sites": state.references["relatedness_sites"],
        "pedigree": [i.model_dump() for i in state.case.pedigree],
        "output_dir": Path.cwd() / "results" / "relatedness",
    }


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
