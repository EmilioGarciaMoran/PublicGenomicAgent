"""PED (validación, parentesco esperado) e IBS local con familias simuladas."""
from __future__ import annotations

import random

import pytest

from publicgenomicagent.tools.pedigree import (
    inbreeding,
    kinship_function,
    parse_ped_text,
    relationship,
    validate_pedigree,
    write_ped,
)
from publicgenomicagent.tools.relatedness import evaluate, pair_stats, parse_genotypes

TRIO = """FAM F 0 0 1 1
FAM M 0 0 2 1
FAM P F M 1 2
"""

# Abuelos comunes -> hermanos A y B -> sus hijos FA y MO (primos) se casan
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


# --- PED -----------------------------------------------------------------

def test_parse_ped_codes_and_unknown_phenotype():
    fids, inds = parse_ped_text("# cabecera\nF1 X 0 0 0 0\n" + TRIO)
    by = {i.sample: i for i in inds}
    assert by["P"].affected and by["P"].father == "F" and by["M"].sex == "F"
    assert by["X"].sex == "U" and not by["X"].phenotype_known
    assert fids == ["F1", "FAM"]


def test_write_ped_roundtrip():
    _, inds = parse_ped_text(CONSANG)
    _, again = parse_ped_text(write_ped(inds, "FAM"))
    assert [i.model_dump() for i in again] == [i.model_dump() for i in inds]


@pytest.mark.parametrize("ped,needle", [
    ("F A 0 0 1 1\nF A 0 0 1 1\n", "duplicados"),
    ("F P X M 1 2\nF M 0 0 2 1\n", "no tiene fila"),
    ("F P F M 1 2\nF F 0 0 2 1\nF M 0 0 2 1\n", "sexo es F"),
    ("F A B 0 1 1\nF B A 0 1 1\n", "ciclo"),
    ("F P P 0 1 2\n", "propio padre"),
])
def test_validation_issues(ped, needle):
    _, inds = parse_ped_text(ped)
    v = validate_pedigree(inds)
    assert not v.ok and any(needle in i for i in v.issues), v.issues


def test_validation_against_bams():
    _, inds = parse_ped_text(TRIO)
    v = validate_pedigree(inds, sequenced={"F", "M", "P", "Z"})
    assert any("no están en el PED" in i for i in v.issues)
    assert validate_pedigree(inds, sequenced={"F", "M", "P"}).ok


def test_kinship_values_standard_and_consanguineous():
    _, trio = parse_ped_text(TRIO)
    phi = kinship_function(trio)
    assert phi("P", "F") == 0.25 and phi("F", "M") == 0 and phi("P", "P") == 0.5
    _, fam = parse_ped_text(CONSANG)
    phi = kinship_function(fam)
    assert phi("FA", "MO") == 1 / 16                    # primos hermanos
    assert inbreeding(fam, "P1") == 1 / 16              # F de hijos de primos
    assert phi("P1", "P2") == pytest.approx(0.25 + 1 / 32)
    assert phi("P1", "FA") == pytest.approx(0.28125)
    assert relationship(fam, "P1", "P2") == "hermanos completos"
    assert relationship(fam, "P1", "FA") == "padre/madre-hijo"
    assert relationship(fam, "FA", "MO") == "otro parentesco"
    assert relationship(fam, "SA", "SB") == "no emparentados"


# --- simulación ----------------------------------------------------------

def simulate(inds, n_sites=2000, seed=7, error=0.0, af=(0.05, 0.5)):
    """Genotipos por transmisión mendeliana; devuelve {sitio: {sample: dosis}}."""
    rng = random.Random(seed)
    by = {i.sample: i for i in inds}
    gts = {}
    for k in range(n_sites):
        p = rng.uniform(*af)
        hap = {}

        def haps(s):
            if s not in hap:
                i = by[s]
                if i.father is None:
                    hap[s] = (int(rng.random() < p), int(rng.random() < p))
                else:
                    hap[s] = (rng.choice(haps(i.father)), rng.choice(haps(i.mother)))
            return hap[s]

        row = {}
        for s in by:
            d = sum(haps(s))
            if error and rng.random() < error:
                d = rng.choice([x for x in (0, 1, 2) if x != d])
            row[s] = d
        gts[("chr1", k + 1)] = row
    return gts


def run(ped_text, gts, seq, **kw):
    _, inds = parse_ped_text(ped_text)
    return evaluate(inds, gts, seq, **kw)


def test_correct_trio_passes_and_parent_child_has_no_ibs0():
    _, inds = parse_ped_text(TRIO)
    status, pairs, issues, _ = run(TRIO, simulate(inds), ["F", "M", "P"])
    assert status == "passed" and not issues
    po = [p for p in pairs if p["declared"] == "padre/madre-hijo"]
    assert len(po) == 2 and all(p["ibs0"] == 0 for p in po)
    assert all(0.4 < p["relatedness"] < 0.6 for p in po)
    fm = next(p for p in pairs if {p["a"], p["b"]} == {"F", "M"})
    assert abs(fm["relatedness"]) < 0.1


def test_wrong_father_fails_without_inferring_cause():
    _, inds = parse_ped_text(TRIO + "FAM X 0 0 1 1\n")
    gts = simulate(inds)
    for row in gts.values():
        row["F"] = row.pop("X")          # el BAM "del padre" es de otra persona
    status, pairs, issues, _ = run(TRIO, gts, ["F", "M", "P"])
    assert status == "failed"
    assert any("incompatible con padre/madre-hijo" in i and "F" in i for i in issues)
    assert not any("paternidad" in i.lower() for i in issues)


def test_duplicate_sample_fails():
    _, inds = parse_ped_text(TRIO)
    gts = simulate(inds)
    for row in gts.values():
        row["P"] = row["M"]
    status, _, issues, _ = run(TRIO, gts, ["F", "M", "P"])
    assert status == "failed" and any("duplicado" in i for i in issues)


def test_consanguineous_family_has_no_false_failures():
    _, inds = parse_ped_text(CONSANG)
    status, pairs, issues, _ = run(CONSANG, simulate(inds, n_sites=3000),
                                   ["FA", "MO", "P1", "P2", "S1"])
    assert status == "passed", issues
    sib = next(p for p in pairs if {p["a"], p["b"]} == {"P1", "P2"})
    assert sib["expected_relatedness"] == pytest.approx(0.5625)


def test_declared_siblings_who_are_not_fail():
    ped = TRIO + "FAM P2 F M 2 2\n"
    _, inds = parse_ped_text(ped + "FAM U 0 0 2 1\n")
    gts = simulate(inds)
    for row in gts.values():
        row["P2"] = row.pop("U")
    status, _, issues, _ = run(ped, gts, ["P", "P2"])
    assert status == "failed" and any("menor que el declarado" in i for i in issues)


def test_undeclared_relationship_is_a_warning():
    true_ped = "FAM G 0 0 1 1\nFAM H 0 0 2 1\nFAM A G H 1 1\nFAM B G H 2 1\n"
    _, true_inds = parse_ped_text(true_ped)
    gts = simulate(true_inds)
    declared = "FAM A 0 0 1 1\nFAM B 0 0 2 2\n"       # declarados no emparentados
    status, _, issues, warnings = run(declared, gts, ["A", "B"])
    assert status == "warning" and not issues
    assert any("no declarado" in w for w in warnings)


def test_genotyping_noise_does_not_break_parent_child():
    _, inds = parse_ped_text(TRIO)
    status, _, issues, _ = run(TRIO, simulate(inds, error=0.002), ["F", "M", "P"])
    assert status == "passed", issues


def test_panel_with_few_sites_is_not_verifiable_not_passed():
    _, inds = parse_ped_text(TRIO)
    status, pairs, _, _ = run(TRIO, simulate(inds, n_sites=20), ["F", "M", "P"])
    assert status == "not_verifiable"
    assert all(p["status"] == "not_verifiable" for p in pairs)


def test_single_sequenced_individual_is_not_verifiable():
    _, inds = parse_ped_text(TRIO)
    status, pairs, _, warnings = run(TRIO, simulate(inds), ["P"])
    assert status == "not_verifiable" and pairs == []


def test_pair_stats_definitions():
    s = pair_stats([0, 2, 1, 1, None], [2, 2, 1, 0, 1], min_hets=1)
    assert s["n"] == 4 and s["ibs0"] == 1 and s["ibs2"] == 2 and s["shared_hets"] == 1


# --- lectura de genotipos de bcftools ------------------------------------

def test_parse_genotypes_filters_depth_and_foreign_alleles():
    sites = [("chr1", 10, "A", "G"), ("chr1", 20, "C", "T"), ("chr1", 30, "G", "A")]
    vcf = [
        "##fileformat=VCFv4.2\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n",
        "chr1\t10\t.\tA\tG\t50\t.\t.\tGT:DP\t0/1:20\t1/1:3\n",      # s2 DP<7
        "chr1\t20\t.\tC\t.\t0\t.\t.\tGT:DP\t0/0:15\t./.:0\n",       # 0/0 llamado
        "chr1\t30\t.\tG\tC\t50\t.\t.\tGT:DP\t0/1:30\t0/0:30\n",     # alelo ajeno
    ]
    samples, g = parse_genotypes(vcf, sites, min_dp=7)
    assert samples == ["s1", "s2"]
    assert g[("chr1", 10)] == {"s1": 1, "s2": None}
    assert g[("chr1", 20)] == {"s1": 0, "s2": None}
    assert g[("chr1", 30)] == {"s1": None, "s2": 0}


# --- CLI: --ped / --relatedness-sites / excepción humana -----------------

def test_cli_family_options(tmp_path):
    import typer
    from publicgenomicagent.agent.state import AgentState
    from publicgenomicagent.cli.main import _apply_family_options
    from publicgenomicagent.tools.base import CaseManifest

    def state():
        st = AgentState(case=CaseManifest(case_id="C", source="clinical", proband="proband",
                                          hpo_terms={"proband": ["HP:0001939"]}))
        st.bams = {"F": tmp_path / "f.bam", "M": tmp_path / "m.bam", "P": tmp_path / "p.bam"}
        return st

    good = tmp_path / "ok.ped"; good.write_text(TRIO)
    sites = tmp_path / "sites.tsv"; sites.write_text("chr1\t10\tA\tG\n")
    st = state()
    _apply_family_options(st, str(good), str(sites), True)
    assert [i.sample for i in st.case.pedigree] == ["F", "M", "P"]
    assert st.case.proband == "P"                         # único afectado
    assert st.references["relatedness_sites"] == sites.resolve()
    assert st.gate_overrides == ["relatedness_not_verifiable"]

    bad = tmp_path / "bad.ped"; bad.write_text("FAM P F M 1 2\nFAM F 0 0 2 1\nFAM M 0 0 2 1\n")
    with pytest.raises(typer.Exit):
        _apply_family_options(state(), str(bad), None, False)
