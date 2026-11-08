"""E1: frecuencia gnomAD v4 y filtro de variantes raras (sin red: gnomAD simulado)."""
from __future__ import annotations

import json

import pytest

import publicgenomicagent.knowledge.gnomad as gn
from publicgenomicagent.tools.frequency import HEADER_MARK, annotate_lines
from publicgenomicagent.tools.mendelian import rare_filter


def gvar(pos, ref="A", alt="G", exome=None, genome=None):
    return {"variant_id": f"1-{pos}-{ref}-{alt}", "pos": pos, "ref": ref, "alt": alt,
            "exome": exome, "genome": genome}


def block(ac, an, mid=(0, 0), totals=True, extra=()):
    pops = [{"id": "mid", "ac": mid[0], "an": mid[1]}, *extra]
    b = {"populations": pops}
    if totals:
        b.update(ac=ac, an=an)
    return b


# --- frecuencias -----------------------------------------------------------

def test_frequencies_sum_exome_and_genome_with_totals():
    v = gvar(10, exome=block(10, 1000, mid=(5, 100)), genome=block(2, 1000, mid=(1, 50)))
    f = gn.variant_frequencies(v)
    assert f["af"] == pytest.approx(12 / 2000) and f["af_mid"] == pytest.approx(6 / 150)
    assert f["method"] == "total"


def test_frequencies_fallback_sums_ancestry_groups_without_sex_splits():
    pops = ({"id": "nfe", "ac": 4, "an": 400}, {"id": "afr", "ac": 1, "an": 500},
            {"id": "XX", "ac": 99, "an": 999}, {"id": "nfe_XX", "ac": 99, "an": 999})
    v = gvar(10, exome=block(0, 0, mid=(2, 100), totals=False, extra=pops))
    f = gn.variant_frequencies(v)
    assert (f["ac"], f["an"]) == (2 + 4 + 1, 100 + 400 + 500)      # mid + nfe + afr
    assert f["method"] == "suma_poblaciones"


def test_frequencies_without_data():
    f = gn.variant_frequencies(gvar(10))
    assert f["af"] is None and f["an"] == 0 and f["method"] == "sin_datos"


# --- anotación del VCF -----------------------------------------------------

HDR = ["##fileformat=VCFv4.2\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\n"]


def vline(pos, ref="A", alt="G", info="."):
    return f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t{info}\tGT\t0/1\n"


def test_annotate_lines_adds_fields_and_header_mark():
    lookup = {("1", 10, "A", "G"): {"af": 0.2, "ac": 20, "an": 100, "af_mid": 0.5, "an_mid": 10}}
    hdr, lines, st = annotate_lines(HDR, [vline(10, info="CLNSIG=x"), vline(20), vline(30, alt="G,T")], lookup)
    assert "CLNSIG=x;PGA_AF=0.2;PGA_AC=20;PGA_AN=100;PGA_AF_MID=0.5;PGA_AN_MID=10" in lines[0]
    assert lines[1] == vline(20)                       # ausente en gnomAD: sin campos
    assert st == {"n": 3, "found": 1, "absent": 1, "multiallelic_skipped": 1, "ref_mismatch": 0}
    assert hdr[-1].startswith("#CHROM") and any(h.startswith(HEADER_MARK) for h in hdr)


# --- filtro de rareza ------------------------------------------------------

ANN_HDR = HDR[:1] + [HEADER_MARK + "gnomad_r4;date=2026-11-08\n"] + HDR[1:]


def test_rare_filter_global_mid_and_absent():
    lines = [vline(1, info="PGA_AF=0.05;PGA_AF_MID=0.06"),        # común global
             vline(2, info="PGA_AF=0.001;PGA_AF_MID=0.03"),       # común solo en Oriente Medio
             vline(3, info="PGA_AF=0.001;PGA_AF_MID=0.002"),      # rara
             vline(4)]                                            # ausente en gnomAD
    kept, st = rare_filter(ANN_HDR, lines, 0.01)
    assert [l.split("\t")[1] for l in kept] == ["3", "4"]
    assert st["removed_common"] == 2 and st["removed_common_mid_only"] == 1
    assert st["kept_absent_in_gnomad"] == 1


def test_rare_filter_refuses_unannotated_vcf():
    with pytest.raises(ValueError, match="no se puede distinguir"):
        rare_filter(HDR, [vline(4)], 0.01)


# --- consulta: ventanas, caché y caída a la consulta básica ----------------

def test_fetch_windows_cache_and_fallback(tmp_path, monkeypatch):
    monkeypatch.setattr(gn, "AF_CACHE", tmp_path)
    calls = []

    def fake_post(query, chrom, a, b, retries=4):
        calls.append((query is gn.AF_QUERY, chrom, a, b))
        if query is gn.AF_QUERY:
            return {"errors": [{"message": "Cannot query field ac"}]}
        return {"data": {"region": {"variants": [gvar(a)]}}}

    monkeypatch.setattr(gn, "_post", fake_post)
    vs = gn.fetch_region_af("chr1", 1, 120_000)
    assert [c[2:] for c in calls if not c[0]] == [(1, 50_000), (50_001, 100_000), (100_001, 120_000)]
    assert all(c[1] == "1" for c in calls) and len(vs) == 3
    calls.clear()
    assert gn.fetch_region_af("1", 1, 120_000) == vs and calls == []      # caché: sin red


# --- orden en el planner ---------------------------------------------------

def test_planner_annotates_frequency_before_mendelian_only_when_requested(tmp_path):
    from publicgenomicagent.agent.decision_planner import DecisionPlanner
    from publicgenomicagent.agent.state import AgentState
    from publicgenomicagent.tools.base import CaseManifest, GenomicRange, IndividualSpec
    from tests.unit.test_decision_planner import _run

    def state(max_af=None, clinvar=False):
        ped = [IndividualSpec(sample="F", sex="M", affected=False),
               IndividualSpec(sample="M", sex="F", affected=False),
               IndividualSpec(sample="P", sex="U", affected=True, father="F", mother="M")]
        st = AgentState(case=CaseManifest(case_id="C", source="clinical", pedigree=ped, proband="P",
                                          candidate_rois=[GenomicRange(chrom="chr1", start=1, end=99)]))
        st.bams = {s: tmp_path / f"{s}.bam" for s in "FMP"}
        st.references["hg38"] = tmp_path / "ref.fa"
        if clinvar:
            st.references["clinvar"] = tmp_path / "cv.vcf.gz"
        st.artifacts["roi_dir"] = tmp_path / "roi"
        if max_af is not None:
            st.params["max_af"] = max_af
        return st

    def drive(st):
        seq, p = [], DecisionPlanner()
        for _ in range(40):
            a = p.next_action(st)
            if a is None:
                return seq, None
            seq.append(a.tool_name)
            if a.tool_name == "mendelian_filter":
                return seq, a
            if a.tool_name == "annotate_frequency":
                st.vcfs["frequency"] = a.args["output_vcf"]
            _run(st, a)
        raise AssertionError(seq)

    seq, mf = drive(state())
    assert "annotate_frequency" not in seq and mf.args["max_af"] is None
    seq, mf = drive(state(0.01, clinvar=True))
    assert seq.index("annotate_variants") < seq.index("annotate_frequency") < seq.index("mendelian_filter")
    assert mf.args["max_af"] == 0.01 and str(mf.args["trio_vcf"]).endswith(".freq.vcf.gz")


def test_mendelian_waits_if_frequency_annotation_never_succeeds(tmp_path):
    from publicgenomicagent.agent.planner import RuleBasedPlanner
    from publicgenomicagent.agent.state import AgentState
    from publicgenomicagent.tools.base import CaseManifest, GenomicRange, IndividualSpec
    ped = [IndividualSpec(sample="F", sex="M", affected=False),
           IndividualSpec(sample="M", sex="F", affected=False),
           IndividualSpec(sample="P", sex="U", affected=True, father="F", mother="M")]
    st = AgentState(case=CaseManifest(case_id="C", source="clinical", pedigree=ped, proband="P",
                                      candidate_rois=[GenomicRange(chrom="chr1", start=1, end=99)]))
    st.params["max_af"] = 0.01
    st.vcfs["trio"] = tmp_path / "joint.vcf.gz"
    assert RuleBasedPlanner()._rule_mendelian(st) is None          # nunca sin frecuencias


def test_ref_mismatch_with_gnomad_is_flagged_not_called_absent():
    """Caso real: Sandbox con referencia sintética (A en 22:42128941; GRCh38 y gnomAD: C)."""
    lookup = {("22", 42128941, "C", "T"): {"af": 0.3, "ac": 3, "an": 10, "af_mid": None, "an_mid": 0},
              ("22", 42128941, "C", "A"): {"af": 0.1, "ac": 1, "an": 10, "af_mid": None, "an_mid": 0}}
    line = "chr22\t42128941\t.\tA\tT\t50\tPASS\t.\tGT\t1/1\n"
    hdr, lines, st = annotate_lines(HDR, [line, vline(99)], lookup)
    assert st["ref_mismatch"] == 1 and st["absent"] == 1
    assert "PGA_REF_MISMATCH" in lines[0] and "PGA_AF" not in lines[0]
    kept, rs = rare_filter(hdr, lines, 0.01)
    assert len(kept) == 2 and rs["kept_ref_mismatch"] == 1 and rs["kept_absent_in_gnomad"] == 1
