"""Segregación recesiva generalizada a familias (hermanos, consanguinidad)."""
from __future__ import annotations

import random

from publicgenomicagent.tools.mendelian import (
    _filter_auto_rec_family,
    _filter_auto_rec_hom,
    family_roles,
)
from publicgenomicagent.tools.pedigree import parse_ped_text

CONSANG = """FAM GF 0 0 1 1
FAM GM 0 0 2 1
FAM A GF GM 1 1
FAM B GF GM 2 1
FAM SA 0 0 2 1
FAM SB 0 0 1 1
FAM FA A SA 1 1
FAM MO SB B 2 1
FAM P1 FA MO 1 2
FAM P2 FA MO 2 2
FAM S1 FA MO 1 1
"""


def line(samples, gts, qual=50, dp=30):
    cols = ["chr1", "100", ".", "A", "G", str(qual), "PASS", ".", "GT:DP"]
    return "\t".join(cols + [f"{gts[s]}:{dp}" for s in samples]) + "\n"


def rec(ped_text, samples, proband, gts, **kw):
    _, ped = parse_ped_text(ped_text)
    roles = family_roles(ped, samples, proband)
    hits, n = _filter_auto_rec_family([line(samples, gts, **kw)], samples, roles, 10, 20)
    return n == 1, roles


SEQ = ["FA", "MO", "P1", "P2", "S1"]          # abuelos y tíos sin BAM


def test_roles_in_consanguineous_family():
    _, ped = parse_ped_text(CONSANG)
    r = family_roles(ped, SEQ, "P1")
    assert r == {"affected": ["P1", "P2"], "carriers": ["FA", "MO"],
                 "carriers_any": [], "unaffected": ["S1"]}


def test_family_rule_consanguineous_with_affected_and_unaffected_sibs():
    ok = {"FA": "0/1", "MO": "0/1", "P1": "1/1", "P2": "1/1", "S1": "0/1"}
    assert rec(CONSANG, SEQ, "P1", ok)[0]
    assert not rec(CONSANG, SEQ, "P1", {**ok, "S1": "1/1"})[0]      # hermano sano homocigoto
    assert not rec(CONSANG, SEQ, "P1", {**ok, "P2": "0/1"})[0]      # afectado no homocigoto
    assert not rec(CONSANG, SEQ, "P1", {**ok, "FA": "0/0"})[0]      # portador obligado sin la variante
    assert rec(CONSANG, SEQ, "P1", {**ok, "S1": "0/0"})[0]
    assert rec(CONSANG, SEQ, "P1", {**ok, "S1": "./."})[0]          # sin dato no descarta
    assert not rec(CONSANG, SEQ, "P1", ok, dp=5)[0]                 # DP bajo en afectados
    assert not rec(CONSANG, SEQ, "P1", ok, qual=10)[0]


def test_siblings_only_without_parents_sequenced():
    seq = ["P1", "P2", "S1"]
    gts = {"P1": "1/1", "P2": "1/1", "S1": "0/1"}
    hit, roles = rec(CONSANG, seq, "P1", gts)
    assert hit and roles["carriers"] == []
    assert not rec(CONSANG, seq, "P1", {**gts, "S1": "1/1"})[0]


def test_affected_parent_pseudodominance_must_be_homozygous():
    ped = "F D 0 0 1 2\nF M 0 0 2 1\nF P D M 1 2\n"           # padre afectado
    seq = ["D", "M", "P"]
    _, inds = parse_ped_text(ped)
    r = family_roles(inds, seq, "P")
    assert r["affected"] == ["P", "D"] and r["carriers"] == ["M"]
    assert rec(ped, seq, "P", {"D": "1/1", "M": "0/1", "P": "1/1"})[0]
    assert not rec(ped, seq, "P", {"D": "0/1", "M": "0/1", "P": "1/1"})[0]


def test_parent_with_unknown_phenotype_may_be_het_or_hom():
    ped = "F D 0 0 1 0\nF M 0 0 2 1\nF P D M 1 2\n"
    seq = ["D", "M", "P"]
    assert family_roles(parse_ped_text(ped)[1], seq, "P")["carriers_any"] == ["D"]
    for d in ("0/1", "1/1"):
        assert rec(ped, seq, "P", {"D": d, "M": "0/1", "P": "1/1"})[0]
    assert not rec(ped, seq, "P", {"D": "0/0", "M": "0/1", "P": "1/1"})[0]


def test_identical_to_classic_trio_rule_on_random_variants():
    """En un trío, la regla familiar da exactamente lo mismo que la clásica.

    Incluye el caso de un probando no marcado como afectado en el PED.
    """
    rng = random.Random(3)
    gts_pool = ["0/0", "0/1", "1/1", "./.", "0|1", "1|1"]
    samples = ["father", "mother", "proband"]
    for affected_flag in ("2", "1"):
        ped = f"F father 0 0 1 1\nF mother 0 0 2 1\nF proband father mother 1 {affected_flag}\n"
        _, inds = parse_ped_text(ped)
        roles = family_roles(inds, samples, "proband")
        lines = []
        for _ in range(2000):
            gts = {s: rng.choice(gts_pool) for s in samples}
            lines.append(line(samples, gts, qual=rng.choice([10, 50]), dp=rng.choice([5, 30])))
        new, _ = _filter_auto_rec_family(lines, samples, roles, 10, 20)
        old, _ = _filter_auto_rec_hom(lines, [], samples, "proband", "father", "mother", 10, 20)
        assert new == old and len(new) > 0


def test_planner_handles_pedigree_members_without_bam(tmp_path):
    """Abuelos y tíos del bucle de consanguinidad no tienen BAM: joint_call debe proponerse."""
    from publicgenomicagent.agent.decision_planner import DecisionPlanner
    from publicgenomicagent.agent.state import AgentState
    from publicgenomicagent.tools.base import CaseManifest, GenomicRange
    from tests.unit.test_decision_planner import _run

    _, ped = parse_ped_text(CONSANG)
    st = AgentState(case=CaseManifest(case_id="F", source="clinical", pedigree=ped, proband="P1",
                                      candidate_rois=[GenomicRange(chrom="chr1", start=1, end=100)]))
    st.bams = {s: tmp_path / f"{s}.bam" for s in SEQ}
    st.references["hg38"] = tmp_path / "ref.fa"
    st.artifacts["roi_dir"] = tmp_path / "roi"
    p = DecisionPlanner()
    seq = []
    for _ in range(60):
        a = p.next_action(st)
        if a is None:
            break
        seq.append(a.tool_name)
        if a.tool_name == "joint_call":
            assert [s for s, _ in a.args["bams"]] == SEQ
        if a.tool_name == "mendelian_filter":
            assert {i["sample"] for i in a.args["pedigree"]} >= {"GF", "GM", "A", "B"}
        _run(st, a)
    assert seq.count("qc_bam") == 5 and "joint_call" in seq and "mendelian_filter" in seq
