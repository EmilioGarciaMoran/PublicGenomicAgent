"""Loop agéntico determinista para análisis de casos genómicos.

Esta primera versión del loop NO usa LLM. Decide la siguiente tool
con reglas explícitas sobre el estado de la sesión. Esto tiene tres
ventajas:

  1. Es testeable hoy, sin dependencia de un modelo externo.
  2. Es auditable: cada decisión tiene una regla legible.
  3. Es el esqueleto sobre el que después se enchufa el LLM:
     el LLM solo tendrá que reemplazar `_decide_next_step`,
     el resto del loop (ejecución, trazabilidad, persistencia)
     se mantiene igual.

El loop se organiza como una secuencia de `Step`s. Cada Step declara:

  - su nombre
  - si sus precondiciones se cumplen (`can_run`)
  - qué tool llamar y con qué argumentos (`plan`)

El `AgentLoop` itera: encuentra el primer step que pueda correr,
lo ejecuta, actualiza estado, y vuelve a empezar. Termina cuando
ningún step puede correr, o cuando se alcanza `max_steps`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from ..tools.base import CaseManifest, GenomicRange
from .planner import Planner
from .state import AgentState, SessionContext


@dataclass
class StepResult:
    """Resultado de intentar ejecutar un step."""

    step_name: str
    tool_name: str | None
    output: Any | None = None
    skipped_reason: str | None = None

    @property
    def executed(self) -> bool:
        return self.tool_name is not None and self.skipped_reason is None


@dataclass
class Step:
    """Un paso del loop: precondición + plan + ejecución.

    - `can_run(state)` decide si el step es aplicable ahora.
    - `plan(state)` devuelve `(tool_name, kwargs)` o `None` si no
      puede planificarse (por ejemplo, faltan rutas que el usuario
      aún no ha proporcionado).
    - `run` se implementa por defecto: llama a `ctx.call_tool` con
      el plan. Los steps que necesiten lógica extra pueden sobrescribir.
    """

    name: str
    can_run: Callable[[AgentState], bool]
    plan: Callable[[AgentState], tuple[str, dict[str, Any]] | None]
    description: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)

    def execute(self, ctx: SessionContext) -> StepResult:
        if not self.can_run(ctx.state):
            return StepResult(
                step_name=self.name,
                tool_name=None,
                skipped_reason=f"precondiciones no cumplidas para '{self.name}'",
            )

        planned = self.plan(ctx.state)
        if planned is None:
            return StepResult(
                step_name=self.name,
                tool_name=None,
                skipped_reason=f"no hay plan ejecutable para '{self.name}'",
            )

        tool_name, kwargs = planned
        out = ctx.call_tool(tool_name, **kwargs)
        return StepResult(
            step_name=self.name,
            tool_name=tool_name,
            output=out,
        )


# ---------------------------------------------------------------------------
# Steps concretos para un caso de trío con ROIs candidatos
# ---------------------------------------------------------------------------

def _samples_in_case(state: AgentState) -> list[str]:
    """Samples del pedigrí del caso (o del dict de BAMs si no hay pedigrí)."""
    if state.case.pedigree:
        return [ind.sample for ind in state.case.pedigree]
    return list(state.bams.keys())


def _rois_for_case(state: AgentState) -> list[GenomicRange]:
    """ROIs candidatos. Si no hay, se derivan de candidate_genes (aún no
    implementado; aquí solo devolvemos lo que haya explícito)."""
    return list(state.case.candidate_rois)


def _has_rois(state: AgentState) -> bool:
    return len(_rois_for_case(state)) > 0


def _has_bams(state: AgentState) -> bool:
    return len(state.bams) > 0


def _has_qc_done(state: AgentState) -> bool:
    samples = set(_samples_in_case(state))
    done = {c.input.get("bam_path") for c in state.calls_of("qc_bam")}
    expected = {str(state.bams[s]) for s in samples if s in state.bams}
    return expected.issubset(done)


def _has_roi_bams(state: AgentState) -> bool:
    """¿Tenemos ya un sub-BAM por ROI por sample?"""
    samples = _samples_in_case(state)
    rois = _rois_for_case(state)
    if not samples or not rois:
        return False
    # Marcamos como hecho si existe un fetch_roi por (sample, roi)
    fetch_calls = state.calls_of("fetch_roi")
    done_pairs = set()
    for c in fetch_calls:
        done_pairs.add((c.input.get("bam_path"), c.input.get("region")))
    for s in samples:
        bam = state.bams.get(s)
        if bam is None:
            return False
        for roi in rois:
            region = f"{roi.chrom}:{roi.start}-{roi.end}"
            if (str(bam), region) not in done_pairs:
                return False
    return True


def _roi_region(roi: GenomicRange) -> str:
    return f"{roi.chrom}:{roi.start}-{roi.end}"


# --- Step: QC de cada BAM ------------------------------------------------

def _plan_qc(state: AgentState) -> tuple[str, dict[str, Any]] | None:
    for sample in _samples_in_case(state):
        bam = state.bams.get(sample)
        if bam is None:
            continue
        already = any(
            c.input.get("bam_path") == str(bam) for c in state.calls_of("qc_bam")
        )
        if not already:
            return ("qc_bam", {"bam_path": bam, "level": "structural"})
    return None


STEP_QC = Step(
    name="qc_bams",
    can_run=lambda s: _has_bams(s) and not _has_qc_done(s),
    plan=_plan_qc,
    description="QC estructural de cada BAM del caso (header, índice, orden).",
    tags=("hts", "qc"),
)


# --- Step: fetch_roi por (sample, roi) ----------------------------------

def _plan_fetch_roi(state: AgentState) -> tuple[str, dict[str, Any]] | None:
    rois = _rois_for_case(state)
    samples = _samples_in_case(state)
    if not rois or not samples:
        return None

    # Recorremos en orden deterministico (sample, roi) y devolvemos
    # el primer par que aún no tenga sub-BAM generado.
    fetch_calls = state.calls_of("fetch_roi")
    done_pairs = {
        (c.input.get("bam_path"), c.input.get("region")) for c in fetch_calls
    }
    for sample in samples:
        bam = state.bams.get(sample)
        if bam is None:
            continue
        for roi in rois:
            region = _roi_region(roi)
            if (str(bam), region) in done_pairs:
                continue
            out_dir = state.artifacts.get("roi_dir") or Path.cwd() / "results" / "roi"
            out_dir = Path(out_dir)
            out_bam = out_dir / f"{sample}.{roi.chrom}_{roi.start}_{roi.end}.bam"
            return (
                "fetch_roi",
                {
                    "bam_path": bam,
                    "region": region,
                    "output_bam": out_bam,
                },
            )
    return None


STEP_FETCH_ROI = Step(
    name="fetch_roi_per_sample_per_roi",
    can_run=lambda s: _has_rois(s) and _has_bams(s) and not _has_roi_bams(s),
    plan=_plan_fetch_roi,
    description="Corta un sub-BAM por (sample, ROI) para variant calling focalizado.",
    tags=("hts", "roi"),
)


# --- Step: call_variants joint por ROI -----------------------------------

def _plan_call_variants(state: AgentState) -> tuple[str, dict[str, Any]] | None:
    rois = _rois_for_case(state)
    samples = _samples_in_case(state)
    if not rois or not samples:
        return None

    # Necesitamos un sub-BAM por sample y ROI, y aún no haber llamado
    # variantes para esa ROI.
    if not _has_roi_bams(state):
        return None

    called_regions = {
        c.input.get("region") for c in state.calls_of("call_variants")
    }

    for roi in rois:
        region = _roi_region(roi)
        if region in called_regions:
            continue

        # Tomamos el primer sample con sub-BAM disponible para esta ROI
        # como "representante" del call (call_variants opera sobre un BAM).
        # En un caso real, el joint calling requiere multi-sample BAM, que
        # es un paso aparte. Aquí invocamos call_variants por ROI sobre el
        # primer sample, y dejamos constancia en las notas.
        for sample in samples:
            bam = state.bams.get(sample)
            if bam is None:
                continue
            out_dir = state.artifacts.get("roi_dir") or Path.cwd() / "results" / "roi"
            sub_bam = Path(out_dir) / f"{sample}.{roi.chrom}_{roi.start}_{roi.end}.bam"
            if not sub_bam.exists():
                continue
            ref = state.references.get("hg38") or state.references.get("ref")
            if ref is None:
                return None  # sin referencia, no se puede llamar
            out_vcf = sub_bam.with_suffix(".vcf.gz")
            return (
                "call_variants",
                {
                    "bam_path": sub_bam,
                    "reference_fasta": ref,
                    "output_vcf": out_vcf,
                    "region": region,
                },
            )
    return None


STEP_CALL_VARIANTS = Step(
    name="call_variants_per_roi",
    can_run=lambda s: _has_roi_bams(s),
    plan=_plan_call_variants,
    description="Llama variantes sobre cada sub-BAM por ROI.",
    tags=("hts", "variants"),
)


# --- Step: mendelian_filter (si hay pedigree y VCF de trío) --------------

def _plan_mendelian(state: AgentState) -> tuple[str, dict[str, Any]] | None:
    if not state.case.has_pedigree():
        return None
    trio_vcf = state.vcfs.get("trio") or state.vcfs.get("joint")
    if trio_vcf is None:
        return None
    if state.has_tool("mendelian_filter"):
        return None  # ya ejecutado

    out_dir = Path.cwd() / "results" / "mendelian"
    return (
        "mendelian_filter",
        {
            "trio_vcf": trio_vcf,
            "output_dir": out_dir,
            "pedigree": [ind.model_dump() for ind in state.case.pedigree],
            "proband": state.case.proband,
        },
    )


STEP_MENDELIAN = Step(
    name="mendelian_filter",
    can_run=lambda s: s.case.has_pedigree(),
    plan=_plan_mendelian,
    description="Filtros mendelianos sobre el VCF del trío.",
    tags=("family", "mendelian"),
)


# ---------------------------------------------------------------------------
# Loop
# ---------------------------------------------------------------------------

DEFAULT_STEPS: list[Step] = [
    STEP_QC,
    STEP_FETCH_ROI,
    STEP_CALL_VARIANTS,
    STEP_MENDELIAN,
]


class AgentLoop:
    """Loop determinista sobre un `SessionContext`.

    Itera pasos hasta que ninguno pueda ejecutarse o se alcance el
    máximo. Deja traza completa en `state.tool_calls` y `state.notes`.
    """

    def __init__(
        self,
        steps: list[Step] | None = None,
        max_steps: int = 50,
        planner: Planner | None = None,
    ):
        self.steps = steps if steps is not None else DEFAULT_STEPS
        self.max_steps = max_steps
        # Si hay planner, decide él. Si no, se usan los steps.
        self.planner = planner

    def run(self, ctx: SessionContext) -> list[StepResult]:
        results: list[StepResult] = []
        for i in range(self.max_steps):
            step_result = self._run_one(ctx)
            if step_result is None:
                ctx.state.note(f"loop terminado tras {i} iteraciones")
                break
            results.append(step_result)
            if step_result.executed and step_result.tool_name:
                ctx.state.note(
                    f"step {step_result.step_name} → {step_result.tool_name} OK"
                )
            elif step_result.skipped_reason:
                ctx.state.note(
                    f"step {step_result.step_name} saltado: "
                    f"{step_result.skipped_reason}"
                )
        else:
            ctx.state.note(f"loop alcanzó max_steps={self.max_steps}")
        return results

    def _run_one(self, ctx: SessionContext) -> StepResult | None:
        """Ejecuta la siguiente acción.

        Si hay `planner`, se le consulta primero. Si devuelve una
        acción, se ejecuta vía `SessionContext.call_tool`. Si el
        planner devuelve None, se cae a los `Step`s clásicos (útil
        para tests y para el modo sin LLM).

        Devuelve `None` si nada aplica.
        """
        if self.planner is not None:
            action = self.planner.next_action(ctx.state)
            if action is not None:
                out = ctx.call_tool(action.tool_name, **action.args)
                return StepResult(
                    step_name=f"planner:{action.tool_name}",
                    tool_name=action.tool_name,
                    output=out,
                )

        for step in self.steps:
            if not step.can_run(ctx.state):
                continue
            planned = step.plan(ctx.state)
            if planned is None:
                continue
            return step.execute(ctx)
        return None


# ---------------------------------------------------------------------------
# Conveniencia: caso canónico de trío
# ---------------------------------------------------------------------------

def build_trio_case(
    *,
    case_id: str,
    father_bam: Path,
    mother_bam: Path,
    proband_bam: Path,
    roi: GenomicRange,
    reference_fasta: Path,
) -> CaseManifest:
    """Construye un `CaseManifest` para un trío clásico.

    Útil para demos y tests. El usuario real normalmente cargará el
    manifiesto desde JSON o dejará que el LLM lo construya.
    """
    from ..tools.base import IndividualSpec

    return CaseManifest(
        case_id=case_id,
        source="clinical",
        pedigree=[
            IndividualSpec(sample="father", sex="M", affected=False),
            IndividualSpec(sample="mother", sex="F", affected=False),
            IndividualSpec(
                sample="proband",
                sex="U",
                affected=True,
                father="father",
                mother="mother",
            ),
        ],
        proband="proband",
        affected_samples=["proband"],
        candidate_rois=[roi],
    )
