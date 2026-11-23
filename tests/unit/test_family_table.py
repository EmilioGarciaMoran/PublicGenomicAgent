"""Tabla familiar para clínicos -> PED, con consanguinidad y aprobación."""
from __future__ import annotations

from fractions import Fraction

import pytest

from publicgenomicagent.tools import family_table as ftm
from publicgenomicagent.tools.family_table import TEMPLATE_EN, TEMPLATE_ES, parse_family_table

HEAD = "id\trelation\tsex\tstatus\tbam\n"


def table(rows: str, parents: str | None = None) -> str:
    return (f"parents_relationship\t{parents}\n" if parents else "") + HEAD + rows


def test_templates_en_es_give_same_pedigree_with_first_cousin_loop():
    en, es = parse_family_table(TEMPLATE_EN), parse_family_table(TEMPLATE_ES)
    assert en.ok and es.ok and en.ped() == es.ped()
    assert en.proband == "P1" and en.expected_f == Fraction(1, 16)
    assert en.computed_f == pytest.approx(1 / 16)
    assert set(en.added) == {"PGF", "PGM", "MGF", "MGM", "FOUNDER1", "FOUNDER2"}


@pytest.mark.parametrize("rel,f", [("unrelated", 0), ("second cousins", 1 / 64),
                                   ("uncle-niece", 1 / 8), ("tío-sobrina", 1 / 8)])
def test_parents_relationship_gives_expected_inbreeding(rel, f):
    ft = parse_family_table(table("P1\tproband\tmale\taffected\t\nFA\tfather\t\tunaffected\t\n"
                                  "MO\tmother\t\tunaffected\t\n", rel))
    assert ft.ok, ft.issues
    assert ft.computed_f == pytest.approx(f)


def test_words_and_spanish_and_deduced_sex():
    ft = parse_family_table(table("Ana\tprobando\tmujer\tafectada\t\nLuis\tpadre\t\tsano\t\n"
                                  "Eva\tmadre\t\tdesconocido\t\nJon\ttío paterno\tvarón\tafectado\t\n"))
    assert ft.ok, ft.issues
    by = {i.sample: i for i in ft.individuals}
    assert by["Luis"].sex == "M" and by["Eva"].phenotype_known is False
    assert by["Jon"].father == by["Luis"].father == "PGF"           # tío = hijo de los abuelos paternos
    assert any("sexo deducido" in n for n in ft.notes)


@pytest.mark.parametrize("rows,msg", [
    ("A\tproband\tmale\taffected\t\nB\tproband\tfemale\taffected\t\n", "exactamente un probando"),
    ("A\tproband\tmale\taffected\t\nB\tneighbour\tmale\taffected\t\n", "relación no reconocida"),
    ("A\tproband\tmale\taffected\t\nB\tfather\tfemale\tunaffected\t\n", "no concuerda con el sexo"),
    ("A\tproband\tmale\tsick-ish\t\n", "estado no reconocido"),
])
def test_clear_errors(rows, msg):
    ft = parse_family_table(table(rows))
    assert not ft.ok and any(msg in e for e in ft.issues)


def test_approval_requires_review_and_detects_changes(tmp_path):
    t = tmp_path / "fam.tsv"
    t.write_text(TEMPLATE_EN)
    with pytest.raises(ValueError, match="revise primero"):
        ftm.approve(t, by="Dr X")
    ftm.record_review(t, drawing=False)
    rec = ftm.approve(t, by="Dr X")
    assert rec["approved_by"] == "Dr X" and rec["expected_inbreeding_proband"] == pytest.approx(1 / 16)
    assert ftm.check_approval(t)["table_sha256"] == rec["table_sha256"]
    t.write_text(TEMPLATE_EN.replace("unaffected", "affected", 1))
    with pytest.raises(ValueError, match="ha cambiado"):
        ftm.check_approval(t)
