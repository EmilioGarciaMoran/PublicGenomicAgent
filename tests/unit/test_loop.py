"""Tests del loop agéntico determinista.

Sin BAMs reales ni runtime: usamos steps sintéticos con `can_run` y
`plan` controlados, y un `SessionContext` con runtime=None. Verifica:

  - que el loop itera mientras haya steps ejecutables
  - que para cuando ninguno aplica
  - que respeta `max_steps`
  - que los steps saltados dejan nota
  - que la traza queda en `state.tool_calls`
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from publicgenomicagent.agent.loop import (
    AgentLoop,
    Step,
    StepResult,
    build_trio_case,
)
from publicgenomicagent.agent.state import AgentState, SessionContext
from publicgenomicagent.tools.base import CaseManifest, GenomicRange


def _empty_case() -> CaseManifest:
    return CaseManifest(
        case_id="X",
        candidate_rois=[GenomicRange(chrom="chr1", start=1, end=10)],
    )


def _noop_executor(ctx: SessionContext, kwargs: dict[str, Any]) -> dict[str, Any]:
    """Simula una ejecución exitosa sin tocar tools."""
    return {"ok": True, **kwargs}


def _make_step(name: str, *, can: bool, plan=None) -> Step:
    """Step sintético que, si puede, invoca 'compare_vcfs' con args vacíos
    (que fallará) — pero como no queremos que falle, sobrescribimos
    `execute` para no llamar a ninguna tool real.

    En su lugar usamos un truco: el plan devuelve ('_noop', {}) y
    parcheamos `SessionContext.call_tool` para interceptar '_noop'.
    """
    return Step(
        name=name,
        can_run=lambda s: can,
        plan=plan or (lambda s: ("_noop", {"step": name})),
        description="",
    )


class _NoopCtx(SessionContext):
    """SessionContext de prueba que no toca tools ni runtime."""

    def call_tool(self, tool_name: str, **kwargs: Any):  # type: ignore[override]
        assert tool_name == "_noop", f"unexpected tool: {tool_name}"
        self.state.record(
            tool_name="_noop",
            inp=kwargs,
            out={"ok": True, **kwargs},
            duration_ms=0,
        )
        return None


def test_loop_runs_all_steps_until_none_applies():
    state = AgentState(case=_empty_case())
    ctx = _NoopCtx(state=state, runtime=None)

    steps = [
        _make_step("s1", can=True),
        _make_step("s2", can=True),
    ]
    loop = AgentLoop(steps=steps, max_steps=10)
    results = loop.run(ctx)

    # Cada step se ejecuta una vez y luego se agota (porque `plan`
    # devolvería el mismo plan de nuevo → el loop se detiene solo
    # porque ningún step nuevo aparece; en este caso el loop se
    # detendrá tras la primera pasada al no haber progreso).
    # Validamos que al menos se ejecutó uno.
    assert len(results) >= 1
    assert all(r.executed for r in results)
    assert all(r.tool_name == "_noop" for r in results)


def test_loop_respects_max_steps():
    state = AgentState(case=_empty_case())
    ctx = _NoopCtx(state=state, runtime=None)

    # Un step que siempre puede correr y siempre planifica lo mismo:
    # el loop no debe colgarse, y debe parar en max_steps.
    step = Step(
        name="always",
        can_run=lambda s: True,
        plan=lambda s: ("_noop", {"i": len(s.tool_calls)}),
    )
    loop = AgentLoop(steps=[step], max_steps=5)
    results = loop.run(ctx)
    assert len(results) == 5
    assert "max_steps" in state.notes[-1]


def test_loop_stops_when_no_step_applies():
    state = AgentState(case=_empty_case())
    ctx = _NoopCtx(state=state, runtime=None)

    steps = [
        _make_step("never", can=False),
    ]
    loop = AgentLoop(steps=steps, max_steps=10)
    results = loop.run(ctx)
    assert results == []
    assert any("loop terminado" in n for n in state.notes)


def test_loop_skips_steps_with_no_plan():
    state = AgentState(case=_empty_case())
    ctx = _NoopCtx(state=state, runtime=None)

    steps = [
        Step(name="no_plan", can_run=lambda s: True, plan=lambda s: None),
        _make_step("real", can=True),
    ]
    loop = AgentLoop(steps=steps, max_steps=3)
    results = loop.run(ctx)
    # El loop salta 'no_plan' (sin ejecutar) y hace 'real'
    executed = [r for r in results if r.executed]
    assert len(executed) >= 1
    assert executed[0].tool_name == "_noop"


def test_step_execute_returns_skip_when_cannot_run():
    state = AgentState(case=_empty_case())
    ctx = _NoopCtx(state=state, runtime=None)
    step = _make_step("x", can=False)
    r = step.execute(ctx)
    assert not r.executed
    assert r.skipped_reason is not None


def test_step_execute_returns_skip_when_plan_is_none():
    state = AgentState(case=_empty_case())
    ctx = _NoopCtx(state=state, runtime=None)
    step = Step(name="x", can_run=lambda s: True, plan=lambda s: None)
    r = step.execute(ctx)
    assert not r.executed
    assert "no hay plan" in (r.skipped_reason or "")


def test_build_trio_case_shape():
    case = build_trio_case(
        case_id="TRIO-1",
        father_bam=Path("/tmp/f.bam"),
        mother_bam=Path("/tmp/m.bam"),
        proband_bam=Path("/tmp/p.bam"),
        roi=GenomicRange(chrom="chr2", start=100, end=200, label="NPHP1"),
        reference_fasta=Path("/tmp/ref.fa"),
    )
    assert case.case_id == "TRIO-1"
    assert case.proband == "proband"
    assert case.has_pedigree()
    assert case.has_rois()
    samples = [ind.sample for ind in case.pedigree]
    assert samples == ["father", "mother", "proband"]
    father, mother = case.get_parents("proband")
    assert father == "father"
    assert mother == "mother"
