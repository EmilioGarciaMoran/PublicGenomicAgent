"""Tests del RuleBasedPlanner, LLMPlanner y los clientes fake."""
from __future__ import annotations

from pathlib import Path

import pytest

from publicgenomicagent.agent.llm import (
    FakeLLMClient,
    LLMResponseFormatError,
    RecordingLLMClient,
    parse_llm_response,
)
from publicgenomicagent.agent.planner import (
    LLMPlanner,
    RuleBasedPlanner,
)
from publicgenomicagent.agent.state import AgentState
from publicgenomicagent.tools.base import (
    CaseManifest,
    GenomicRange,
    IndividualSpec,
)


def _case_with_roi() -> CaseManifest:
    return CaseManifest(
        case_id="T1",
        source="clinical",
        candidate_rois=[GenomicRange(chrom="chr1", start=1, end=100, label="TEST")],
    )


def _trio_case() -> CaseManifest:
    return CaseManifest(
        case_id="TRIO",
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
        candidate_rois=[GenomicRange(chrom="chr1", start=1, end=100)],
    )


# ---------------------------------------------------------------------------
# parse_llm_response
# ---------------------------------------------------------------------------

def test_parse_valid_tool_response():
    r = parse_llm_response(
        '{"tool_name": "qc_bam", "args": {"bam_path": "x.bam"}, '
        '"rationale": "QC pendiente"}'
    )
    assert r.tool_name == "qc_bam"
    assert r.args == {"bam_path": "x.bam"}
    assert r.stop is False


def test_parse_stop_response():
    r = parse_llm_response('{"stop": true, "rationale": "nada que hacer"}')
    assert r.stop is True
    assert r.tool_name is None


def test_parse_tolerates_json_fence():
    r = parse_llm_response(
        '```json\n{"tool_name": "qc_bam", "args": {}, "rationale": ""}\n```'
    )
    assert r.tool_name == "qc_bam"


def test_parse_rejects_empty():
    with pytest.raises(LLMResponseFormatError):
        parse_llm_response("")


def test_parse_rejects_non_json():
    with pytest.raises(LLMResponseFormatError):
        parse_llm_response("not json at all")


def test_parse_rejects_missing_tool_name_and_stop():
    with pytest.raises(LLMResponseFormatError):
        parse_llm_response('{"args": {}}')


def test_parse_rejects_non_object_args():
    with pytest.raises(LLMResponseFormatError):
        parse_llm_response(
            '{"tool_name": "qc_bam", "args": [1,2,3], "rationale": ""}'
        )


# ---------------------------------------------------------------------------
# RuleBasedPlanner
# ---------------------------------------------------------------------------

def test_rule_based_returns_qc_first():
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/tmp/p.bam")}

    p = RuleBasedPlanner()
    action = p.next_action(state)
    assert action is not None
    assert action.tool_name == "qc_bam"
    assert "proband" in action.rationale


def test_rule_based_returns_qc_when_no_rois():
    state = AgentState(case=CaseManifest(case_id="x"))
    state.bams = {"proband": Path("/tmp/p.bam")}

    p = RuleBasedPlanner()
    action = p.next_action(state)
    # Solo hay QC → la primera acción será qc_bam, no fetch_roi.
    assert action is not None
    assert action.tool_name == "qc_bam"


def test_rule_based_skips_qc_if_already_done():
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/tmp/p.bam")}
    state.record(
        tool_name="qc_bam",
        inp={"bam_path": "/tmp/p.bam", "level": "structural"},
        out={"ok": True},
        duration_ms=1,
    )

    p = RuleBasedPlanner()
    action = p.next_action(state)
    assert action is not None
    assert action.tool_name == "fetch_roi"


def test_rule_based_returns_none_when_everything_done():
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/tmp/p.bam")}
    state.record(
        tool_name="qc_bam",
        inp={"bam_path": "/tmp/p.bam"},
        out={"ok": True},
        duration_ms=1,
    )
    state.record(
        tool_name="fetch_roi",
        inp={
            "bam_path": "/tmp/p.bam",
            "region": "chr1:1-100",
            "output_bam": "/tmp/roi.bam",
        },
        out={"ok": True},
        duration_ms=1,
    )
    # Sin referencia, call_variants no aplica; sin VCF de trío,
    # mendelian_filter tampoco. Debe devolver None.
    state.references = {}

    p = RuleBasedPlanner()
    action = p.next_action(state)
    assert action is None


def test_rule_based_trio_returns_mendelian_when_vcf_present():
    state = AgentState(case=_trio_case())
    state.bams = {
        "father": Path("/tmp/f.bam"),
        "mother": Path("/tmp/m.bam"),
        "proband": Path("/tmp/p.bam"),
    }
    state.record(
        tool_name="mendelian_filter",
        inp={},
        out={"ok": True},
        duration_ms=1,
    )
    # Ya ejecutado → no debe proponerlo otra vez
    p = RuleBasedPlanner()
    action = p.next_action(state)
    # Puede proponer QC o fetch_roi antes, pero nunca mendelian_filter
    if action is not None:
        assert action.tool_name != "mendelian_filter"


# ---------------------------------------------------------------------------
# LLMPlanner — con FakeLLMClient
# ---------------------------------------------------------------------------

def test_llm_planner_uses_client_response_when_valid():
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/tmp/p.bam")}

    client = FakeLLMClient(
        [
            '{"tool_name": "qc_bam", '
            '"args": {"bam_path": "/tmp/p.bam", "level": "structural"}, '
            '"rationale": "QC pendiente"}'
        ]
    )
    planner = LLMPlanner(client=client)
    action = planner.next_action(state)

    assert action is not None
    assert action.tool_name == "qc_bam"
    assert action.args["bam_path"] == "/tmp/p.bam"
    assert len(client.calls) == 1


def test_llm_planner_returns_none_on_stop():
    state = AgentState(case=_case_with_roi())
    client = FakeLLMClient(['{"stop": true, "rationale": "nada"}'])
    planner = LLMPlanner(client=client)
    action = planner.next_action(state)
    assert action is None
    assert any("stop" in n for n in state.notes)


def test_llm_planner_falls_back_on_unknown_tool():
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/tmp/p.bam")}

    client = FakeLLMClient(
        ['{"tool_name": "hack_the_planet", "args": {}, "rationale": "..."}']
    )
    planner = LLMPlanner(client=client)
    action = planner.next_action(state)

    # Debe caer al fallback determinista → qc_bam
    assert action is not None
    assert action.tool_name == "qc_bam"
    assert any("tool desconocida" in n for n in state.notes)


def test_llm_planner_falls_back_on_bad_json():
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/tmp/p.bam")}

    client = FakeLLMClient(["esto no es json"])
    planner = LLMPlanner(client=client)
    action = planner.next_action(state)

    assert action is not None
    assert action.tool_name == "qc_bam"
    assert any("formato inválido" in n for n in state.notes)


def test_llm_planner_falls_back_when_client_raises():
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/tmp/p.bam")}

    client = FakeLLMClient([])  # agotado desde el primer call
    planner = LLMPlanner(client=client)
    action = planner.next_action(state)

    assert action is not None
    assert action.tool_name == "qc_bam"
    assert any("error del cliente" in n for n in state.notes)


def test_recording_client_captures_prompts():
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/tmp/p.bam")}

    client = RecordingLLMClient(
        '{"tool_name": "qc_bam", '
        '"args": {"bam_path": "/tmp/p.bam", "level": "structural"}, '
        '"rationale": ""}'
    )
    planner = LLMPlanner(client=client)
    planner.next_action(state)

    assert len(client.calls) == 1
    assert "ESTADO ACTUAL" in client.calls[0]["user"]
    assert "CATÁLOGO DE TOOLS" in client.calls[0]["user"]


def test_llm_planner_does_not_send_paths_by_default():
    """Privacidad: por defecto no se envían rutas absolutas al LLM."""
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/home/hospital/paciente_1234.exoma.bam")}

    client = RecordingLLMClient(
        '{"tool_name": "qc_bam", '
        '"args": {"bam_path": "/home/hospital/paciente_1234.exoma.bam", '
        '"level": "structural"}, '
        '"rationale": ""}'
    )
    planner = LLMPlanner(client=client)
    planner.next_action(state)

    user_prompt = client.calls[0]["user"]
    # El directorio absoluto NO debe aparecer en el prompt.
    assert "/home/hospital/" not in user_prompt
    # Pero el basename sí, para que el LLM sepa qué ficheros hay.
    assert "paciente_1234.exoma.bam" in user_prompt


def test_llm_planner_sends_paths_when_explicitly_enabled():
    state = AgentState(case=_case_with_roi())
    state.bams = {"proband": Path("/home/hospital/paciente_1234.exoma.bam")}

    client = RecordingLLMClient('{"stop": true, "rationale": ""}')
    planner = LLMPlanner(client=client, include_paths=True)
    planner.next_action(state)

    user_prompt = client.calls[0]["user"]
    assert "/home/hospital/" in user_prompt
