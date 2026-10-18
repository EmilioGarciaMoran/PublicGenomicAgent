"""Tests unitarios de las funciones de post-proceso de NL→manifest.

Los tests end-to-end con LLM mockeado no son viables sin refactorizar
la firma de nl_to_case_manifest (dispatch pasa runtime, inp y no
permite un tercer argumento). En su lugar, los tests cubren:

  - _extract_json_block: parseo robusto de respuestas del LLM.
  - _validate_and_clean: post-proceso completo del manifest.

El pipeline end-to-end con LLM real se ha verificado manualmente
(ver docs/backlog.md, bloque A).
"""
from __future__ import annotations

import json
import pytest

from publicgenomicagent.tools.case_manifest import (
    _extract_json_block,
    _sanitize_lists,
    _validate_and_clean,
)


# ---------------------------------------------------------------------------
# _extract_json_block
# ---------------------------------------------------------------------------

def test_extract_json_clean() -> None:
    assert _extract_json_block('{"a": 1}') == {"a": 1}


def test_extract_json_fences() -> None:
    raw = "```json\n{\"a\": 1}\n```"
    assert _extract_json_block(raw) == {"a": 1}


def test_extract_json_embedded() -> None:
    raw = 'Aquí va: {"a": 1} y más texto.'
    assert _extract_json_block(raw) == {"a": 1}


def test_extract_json_invalid_raises() -> None:
    with pytest.raises(ValueError, match="JSON"):
        _extract_json_block('{"a": 1, "b": [')


# ---------------------------------------------------------------------------
# _validate_and_clean: post-proceso completo
# ---------------------------------------------------------------------------

def _base_data() -> dict:
    return {
        "case_id": "t",
        "pedigree": [
            {"sample": "father", "sex": "M", "affected": False,
             "father": None, "mother": None},
            {"sample": "proband", "sex": "M", "affected": True,
             "father": "father", "mother": None},
        ],
        "proband": "proband",
        "affected_samples": ["proband"],
        "consanguinity": False,
    }


def test_normalize_affected_null_to_false() -> None:
    data = {
        "case_id": "t",
        "pedigree": [
            {"sample": "father", "sex": "M", "affected": None,
             "father": None, "mother": None},
        ],
        "proband": None,
        "affected_samples": [],
        "consanguinity": False,
    }
    manifest, _ = _validate_and_clean(data, "text")
    assert manifest.pedigree[0].affected is False


def test_sanitize_removes_null_samples() -> None:
    data = {
        "case_id": "t",
        "pedigree": [
            {"sample": "null", "sex": "U", "affected": False,
             "father": None, "mother": None},
            {"sample": "proband", "sex": "M", "affected": True,
             "father": None, "mother": None},
        ],
        "proband": "proband",
        "affected_samples": ["proband"],
        "consanguinity": False,
    }
    manifest, _ = _validate_and_clean(data, "text")
    samples = [ind.sample for ind in manifest.pedigree]
    assert "null" not in samples
    assert "proband" in samples


def test_sanitize_removes_self_reference() -> None:
    data = {
        "case_id": "t",
        "pedigree": [
            {"sample": "father", "sex": "M", "affected": False,
             "father": "father", "mother": None},
        ],
        "proband": None,
        "affected_samples": [],
        "consanguinity": False,
    }
    manifest, _ = _validate_and_clean(data, "text")
    assert manifest.pedigree[0].father is None


def test_affected_samples_derived_from_pedigree() -> None:
    data = _base_data()
    data["affected_samples"] = ["proband", "father"]  # LLM mete a father por error
    manifest, warnings = _validate_and_clean(data, "text")
    assert manifest.affected_samples == ["proband"]
    assert any("recalculado" in w for w in warnings)


def test_pedigree_consistency_adds_referenced() -> None:
    data = {
        "case_id": "t",
        "pedigree": [
            {"sample": "proband", "sex": "M", "affected": True,
             "father": "unknown_father", "mother": None},
        ],
        "proband": "proband",
        "affected_samples": ["proband"],
        "consanguinity": False,
    }
    manifest, warnings = _validate_and_clean(data, "text")
    samples = [ind.sample for ind in manifest.pedigree]
    assert "unknown_father" in samples
    assert any("unknown_father" in w for w in warnings)


def test_sanitize_lists_cuts_hpo_hallucination() -> None:
    data = {"hpo_terms": {"proband": [f"HP:{i:07d}" for i in range(50)]}}
    data, warnings = _sanitize_lists(data)
    assert data["hpo_terms"]["proband"] == []
    assert any("hpo_terms" in w for w in warnings)


def test_sanitize_lists_cuts_candidate_genes() -> None:
    data = {"candidate_genes": [f"GENE{i}" for i in range(20)]}
    data, warnings = _sanitize_lists(data)
    assert data["candidate_genes"] == []
    assert any("candidate_genes" in w for w in warnings)
