"""Tests de la puerta de QC: ningún análisis sin verificación aprobada."""
from __future__ import annotations

from pathlib import Path

from publicgenomicagent.agent.decision_planner import DecisionPlanner
from publicgenomicagent.agent.loop import AgentLoop
from publicgenomicagent.agent.planner import PlannedAction, RuleBasedPlanner
from publicgenomicagent.agent.qc_gate import qc_args, qc_gate
from publicgenomicagent.agent.state import AgentState, SessionContext
from publicgenomicagent.tools.base import (
    CaseManifest,
    GenomicRange,
    IndividualSpec,
    QCHeader,
)
from publicgenomicagent.tools.qc import compare_reference


def _trio_state(tmp: Path) -> AgentState:
    case = CaseManifest(
        case_id="G", source="clinical",
        pedigree=[
            IndividualSpec(sample="father", sex="M", affected=False),
            IndividualSpec(sample="mother", sex="F", affected=False),
            IndividualSpec(sample="proband", sex="U", affected=True,
                           father="father", mother="mother"),
        ],
        proband="proband",
        candidate_rois=[GenomicRange(chrom="chr1", start=1, end=100)],
    )
    st = AgentState(case=case)
    st.bams = {s: tmp / f"{s}.bam" for s in ("father", "mother", "proband")}
    st.references["hg38"] = tmp / "ref.fa"
    st.artifacts["roi_dir"] = tmp / "roi"
    return st


def _qc(st: AgentState, sample: str, passed: bool | None, ok: bool = True) -> None:
    out = None if passed is None else {
        "passed": passed, "issues": [] if passed else [f"SM esperado '{sample}' no está"]}
    st.record(tool_name="qc_bam", inp={"bam_path": str(st.bams[sample])},
              out=out, duration_ms=1, error=None if ok else "boom")


# --- estados de la puerta ------------------------------------------------

def test_gate_pending_until_every_bam_has_qc(tmp_path):
    st = _trio_state(tmp_path)
    _qc(st, "father", True)
    g = qc_gate(st)
    assert not g.open and set(g.pending) == {"mother", "proband"}


def test_gate_open_when_all_passed(tmp_path):
    st = _trio_state(tmp_path)
    for s in ("father", "mother", "proband"):
        _qc(st, s, True)
    assert qc_gate(st).open


def test_gate_closed_on_failed_qc_with_reason(tmp_path):
    st = _trio_state(tmp_path)
    _qc(st, "father", True); _qc(st, "mother", False); _qc(st, "proband", True)
    g = qc_gate(st)
    assert not g.open and "mother" in g.failed
    assert "mother" in g.reason()


def test_gate_closed_when_output_lacks_passed(tmp_path):
    st = _trio_state(tmp_path)
    for s in ("father", "mother", "proband"):
        st.record(tool_name="qc_bam", inp={"bam_path": str(st.bams[s])},
                  out={"ok": True}, duration_ms=1)
    assert not qc_gate(st).open


def test_crashed_qc_counts_as_pending(tmp_path):
    st = _trio_state(tmp_path)
    _qc(st, "father", None, ok=False)
    assert "father" in qc_gate(st).pending


def test_case_without_bams_has_nothing_to_verify(tmp_path):
    st = AgentState(case=CaseManifest(case_id="P", source="clinical",
                                      hpo_terms={"x": ["HP:0001939"]}))
    assert qc_gate(st).open


def test_qc_args_request_counts_sample_and_reference(tmp_path):
    st = _trio_state(tmp_path)
    (tmp_path / "ref.fa.fai").write_text("chr1\t1000\t6\t60\t61\n")
    a = qc_args(st, "mother", st.bams["mother"])
    assert a["level"] == "counts" and a["expected_sample"] == "mother"
    assert a["expected_reference_fai"] == tmp_path / "ref.fa.fai"


# --- los planners no proponen análisis con la puerta cerrada ---------------

def test_rule_planner_stops_after_failed_qc(tmp_path):
    st = _trio_state(tmp_path)
    _qc(st, "father", True); _qc(st, "mother", False); _qc(st, "proband", True)
    assert RuleBasedPlanner().next_action(st) is None
    assert any("puerta QC cerrada" in n for n in st.notes)


def test_rule_planner_runs_qc_for_all_before_analysis(tmp_path):
    st = _trio_state(tmp_path)
    _qc(st, "father", True)
    a = RuleBasedPlanner().next_action(st)
    assert a.tool_name == "qc_bam" and a.args["expected_sample"] == "mother"


def test_decision_planner_stops_after_failed_qc(tmp_path):
    st = _trio_state(tmp_path)
    _qc(st, "father", True); _qc(st, "mother", False); _qc(st, "proband", True)
    assert DecisionPlanner().next_action(st) is None
    assert any("puerta QC cerrada" in n for n in st.notes)


# --- el loop es el punto único de control --------------------------------

class _PushyPlanner:
    """Planner (o LLM) que propone análisis aunque la puerta esté cerrada."""

    def next_action(self, state):
        return PlannedAction(tool_name="fetch_roi", args={
            "bam_path": "/x.bam", "region": "chr1:1-100", "output_bam": "/y.bam"})


def test_loop_refuses_analysis_from_any_planner_when_gate_closed(tmp_path):
    st = _trio_state(tmp_path)
    _qc(st, "father", True); _qc(st, "mother", False); _qc(st, "proband", True)
    before = len(st.tool_calls)
    results = AgentLoop(planner=_PushyPlanner()).run(SessionContext(state=st))
    assert results == [] and len(st.tool_calls) == before
    assert any("loop: puerta QC cerrada; 'fetch_roi'" in n for n in st.notes)


def test_loop_classic_steps_also_blocked(tmp_path):
    st = _trio_state(tmp_path)
    _qc(st, "father", True); _qc(st, "mother", False); _qc(st, "proband", True)
    before = len(st.tool_calls)
    AgentLoop(planner=None).run(SessionContext(state=st))
    assert len(st.tool_calls) == before          # ni fetch_roi ni nada después


# --- referencia: nombres y longitudes -----------------------------------

def test_compare_reference_flags_length_mismatch_as_other_build():
    h = QCHeader(references=["chr1"], reference_lengths={"chr1": 249250621})  # hg19
    issues, _ = compare_reference(h, {"chr1": 248956422})                     # hg38
    assert issues and "otra build" in issues[0]


def test_compare_reference_ok_and_missing_contig():
    h = QCHeader(references=["chr1", "chrUn"], reference_lengths={"chr1": 10, "chrUn": 5})
    issues, _ = compare_reference(h, {"chr1": 10})
    assert len(issues) == 1 and "ausentes" in issues[0]
    assert compare_reference(QCHeader(references=["chr1"], reference_lengths={"chr1": 10}),
                             {"chr1": 10}) == ([], [])
