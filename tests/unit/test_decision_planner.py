"""Tests del DecisionPlanner sobre el AgentState real."""
from __future__ import annotations

import json

from publicgenomicagent.agent.decision_planner import DecisionPlanner, OPTIONAL_TOOLS, SKIP
from publicgenomicagent.agent.planner import RuleFirstLLMPlanner
from publicgenomicagent.agent.state import AgentState, ToolCallRecord
from publicgenomicagent.tools.base import CaseManifest, GenomicRange, IndividualSpec


def FakeState(tmp, trio=True, text=False, clinvar=True):
    ped = [
        IndividualSpec(sample="F", sex="M", affected=False),
        IndividualSpec(sample="M", sex="F", affected=False),
        IndividualSpec(sample="P", sex="U", affected=True, father="F", mother="M"),
    ] if trio else []
    case = CaseManifest(
        case_id="T", source="clinical", pedigree=ped,
        proband="P" if trio else None,
        candidate_rois=[GenomicRange(chrom="chr1", start=100, end=200)],
        phenotype_text={"P": "renal cysts"} if text else {},
    )
    st = AgentState(case=case)
    samples = ("P", "M", "F") if trio else ("P",)
    st.bams = {s: tmp / f"{s}.bam" for s in samples}
    st.references["hg38"] = tmp / "ref.fa"
    if clinvar:
        st.references["clinvar"] = tmp / "clinvar.vcf.gz"
    st.artifacts["roi_dir"] = tmp / "roi"
    return st


def _run(st, a, ok=True):
    """Simula la ejecución de una acción: registra la llamada y sus efectos."""
    inp = {k: (v if isinstance(v, (list, dict)) else str(v)) for k, v in a.args.items()}
    st.tool_calls.append(ToolCallRecord(tool_name=a.tool_name, input=inp, ok=ok))
    if not ok:
        return
    (st.artifacts["roi_dir"]).mkdir(exist_ok=True)
    if a.tool_name == "fetch_roi": a.args["output_bam"].touch()
    if a.tool_name == "joint_call": st.vcfs["trio"] = a.args["output_vcf"]
    if a.tool_name == "annotate_variants": st.vcfs["annotated"] = a.args["output_vcf"]
    if a.tool_name == "prioritize_variants": st.vcfs["prioritized"] = a.args["output_vcf"]
    if a.tool_name == "mendelian_filter": st.outputs["mendelian_filter"] = {}
    if a.tool_name == "extract_hpo": st.outputs["extract_hpo"] = {"hpo_ids": ["HP:0000113"]}
    if a.tool_name == "phenotype_ranking": st.outputs["phenotype_ranking"] = {}
    if a.tool_name == "prioritize_variants": st.outputs["prioritize_variants"] = {}

class Chooser:
    model_id = "fake"
    def __init__(self, reply=None, boom=False): self.reply, self.boom, self.calls = reply, boom, []
    def choose(self, system, user, schema):
        self.calls.append(schema)
        if self.boom: raise TimeoutError("ollama down")
        return self.reply(schema) if callable(self.reply) else self.reply

def drive(p, st, limit=40, ok=True):
    seq = []
    for _ in range(limit):
        a = p.next_action(st)
        if a is None: return seq
        seq.append(a.tool_name); _run(st, a, ok)
    raise AssertionError(f"no termina: {seq}")

def test_trio_pipeline_order_and_termination_without_llm(tmp_path):
    st = FakeState(tmp_path); ch = Chooser(lambda s: json.dumps({"choice": s["properties"]["choice"]["enum"][0]}))
    p = DecisionPlanner(chooser=ch, skippable=OPTIONAL_TOOLS)
    seq = drive(p, st)
    assert seq.count("qc_bam") == 3 and seq.count("fetch_roi") == 3
    assert seq.index("joint_call") < seq.index("annotate_variants") < seq.index("mendelian_filter") < seq.index("prioritize_variants")
    assert len(seq) == 13   # qc*3 + fetch*3 + call*3 + joint + annotate + mendelian + prioritize
    assert seq[-1] == "prioritize_variants"

def test_llm_never_consulted_when_no_candidates(tmp_path):
    st = FakeState(tmp_path); ch = Chooser(lambda s: json.dumps({"choice": SKIP}))
    p = DecisionPlanner(chooser=ch, skippable=OPTIONAL_TOOLS); drive(p, st)
    assert p.next_action(st) is None
    for schema in ch.calls:   # solo se llamó en bifurcaciones, siempre con enum acotado
        assert SKIP in schema["properties"]["choice"]["enum"]

def test_llm_not_consulted_for_required_steps(tmp_path):
    st = FakeState(tmp_path); ch = Chooser(lambda s: json.dumps({"choice": SKIP}))
    p = DecisionPlanner(chooser=ch, skippable=OPTIONAL_TOOLS)
    for _ in range(9): _run(st, p.next_action(st))     # qc*3 + fetch*3 + call*3
    assert ch.calls == []

def test_skip_skips_optionals_but_not_required(tmp_path):
    st = FakeState(tmp_path); ch = Chooser(json.dumps({"choice": SKIP}))
    p = DecisionPlanner(chooser=ch, skippable=OPTIONAL_TOOLS)
    seq = drive(p, st)
    assert any("SKIP annotate_variants" in n for n in st.notes)   # auditable
    assert "annotate_variants" not in seq and "prioritize_variants" not in seq
    assert "mendelian_filter" in seq and "joint_call" in seq

def test_llm_failure_fails_closed_to_default(tmp_path):
    st = FakeState(tmp_path); p = DecisionPlanner(chooser=Chooser(boom=True), skippable=OPTIONAL_TOOLS)
    seq = drive(p, st)
    assert "annotate_variants" in seq and "prioritize_variants" in seq     # = comportamiento determinista
    assert any(d["event"] == "llm_fallback" for d in p.decisions)

def test_out_of_enum_reply_is_rejected(tmp_path):
    st = FakeState(tmp_path); p = DecisionPlanner(chooser=Chooser(json.dumps({"choice": "rm_rf"})), skippable=OPTIONAL_TOOLS)
    assert "annotate_variants" in drive(p, st)
    assert any(d["event"] == "llm_fallback" for d in p.decisions)

def test_garbage_json_is_rejected(tmp_path):
    st = FakeState(tmp_path); p = DecisionPlanner(chooser=Chooser("blah ```json"), skippable=OPTIONAL_TOOLS)
    assert "prioritize_variants" in drive(p, st)

def test_no_chooser_equals_deterministic(tmp_path):
    a = drive(DecisionPlanner(), FakeState(tmp_path / "a" if (tmp_path/"a").mkdir() is None else tmp_path))
    assert a[0] == "qc_bam"

def test_failed_tool_not_retried_forever(tmp_path):
    st = FakeState(tmp_path); p = DecisionPlanner()
    seq = drive(p, st, ok=False)
    assert seq.count("qc_bam") <= 2 * 3         # max_attempts por (tool,args)
    assert len(seq) < 40

def test_hpo_text_adds_optional_branch(tmp_path):
    st = FakeState(tmp_path, text=True)
    seq = drive(DecisionPlanner(), st)
    assert seq.index("extract_hpo") < seq.index("phenotype_ranking")

def test_decisions_are_logged_with_prompt_hash(tmp_path):
    st = FakeState(tmp_path); p = DecisionPlanner(chooser=Chooser(lambda s: json.dumps({"choice": s["properties"]["choice"]["enum"][0]})), skippable=OPTIONAL_TOOLS)
    drive(p, st)
    llm = [d for d in p.decisions if d["event"] == "llm_choice"]
    assert llm and all(d["model"] == "fake" and len(d["prompt_sha256"]) == 16 for d in llm)
    assert any("decision_planner[" in n for n in st.notes)

def test_old_hybrid_planner_consults_llm_at_the_end(tmp_path):
    """Regresión documentada: el planner actual pregunta al LLM cuando ya no hay nada que hacer."""
    class C:
        n = 0
        def complete(self, s, u): C.n += 1; return '{"stop": true, "rationale": "x"}'
    st = FakeState(tmp_path); p = RuleFirstLLMPlanner(client=C())
    for _ in range(40):
        a = p.next_action(st)
        if a is None: break
        _run(st, a)
    assert C.n >= 1


def test_ollama_chooser_sends_enum_schema_and_parses_reply(tmp_path):
    import threading
    from http.server import BaseHTTPRequestHandler, HTTPServer
    from publicgenomicagent.agent.decision_planner import OllamaChooser

    seen = {}

    class H(BaseHTTPRequestHandler):
        def do_POST(self):
            seen["body"] = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            out = json.dumps({"response": json.dumps({"choice": SKIP})}).encode()
            self.send_response(200); self.send_header("Content-Length", str(len(out)))
            self.end_headers(); self.wfile.write(out)
        def do_GET(self):
            out = json.dumps({"models": [{"name": "qwen2.5:3b", "digest": "abcdef1234567890"}]}).encode()
            self.send_response(200); self.send_header("Content-Length", str(len(out)))
            self.end_headers(); self.wfile.write(out)
        def log_message(self, *a): pass

    srv = HTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        ch = OllamaChooser(model="qwen2.5:3b",
                           endpoint=f"http://127.0.0.1:{srv.server_port}/api/generate", timeout=5)
        st = FakeState(tmp_path)
        p = DecisionPlanner(chooser=ch, skippable=OPTIONAL_TOOLS)
        seq = drive(p, st)
    finally:
        srv.shutdown()
    fmt = seen["body"]["format"]
    assert SKIP in fmt["properties"]["choice"]["enum"]
    assert seen["body"]["options"]["temperature"] == 0
    assert "annotate_variants" not in seq            # SKIP respetado
    assert ch.model_id == "qwen2.5:3b@abcdef123456"
    assert any(d.get("model", "").endswith("@abcdef123456") for d in p.decisions)


def test_default_never_asks_llm_and_runs_hpo_step_even_if_llm_says_skip(tmp_path):
    """Regresión: el 3B descartó phenotype_ranking con HPO presente (sesión real CYP2D6)."""
    st = FakeState(tmp_path, clinvar=False)
    st.case.hpo_terms = {"P": ["HP:0001939"]}
    ch = Chooser(json.dumps({"choice": SKIP}))
    seq = drive(DecisionPlanner(chooser=ch), st)
    assert "phenotype_ranking" in seq and "prioritize_variants" in seq
    assert seq.index("phenotype_ranking") < seq.index("prioritize_variants")
    assert ch.calls == []


def test_only_declared_skippable_steps_are_offered(tmp_path):
    st = FakeState(tmp_path, clinvar=False)
    st.case.hpo_terms = {"P": ["HP:0001939"]}
    ch = Chooser(json.dumps({"choice": SKIP}))
    p = DecisionPlanner(chooser=ch, skippable=frozenset({"prioritize_variants"}))
    seq = drive(p, st)
    assert "phenotype_ranking" in seq and "prioritize_variants" not in seq
    assert len(ch.calls) == 1
    assert ch.calls[0]["properties"]["choice"]["enum"] == ["prioritize_variants", SKIP]
    assert any("SKIP prioritize_variants" in n for n in st.notes)


def test_parse_skippable_rejects_required_steps():
    from publicgenomicagent.agent.decision_planner import parse_skippable
    import pytest
    assert parse_skippable("") == frozenset()
    assert parse_skippable("prioritize_variants, annotate_variants") == frozenset(
        {"prioritize_variants", "annotate_variants"})
    with pytest.raises(ValueError):
        parse_skippable("mendelian_filter")
