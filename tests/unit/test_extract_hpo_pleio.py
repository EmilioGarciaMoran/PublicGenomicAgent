"""Unit test for the pleio-hpo backend of extract_hpo.

Marked as integration because it requires the pga-hpo-cpu
environment and the pleio-hpo assets (~570 MB). Skips gracefully
if the environment is not available.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from publicgenomicagent.env.micromamba import bin_path
from publicgenomicagent.tools.base import ExtractHPOInput
from publicgenomicagent.tools.extract_hpo import extract_hpo


pytestmark = pytest.mark.integration


def _pleio_env_available() -> bool:
    return bin_path("pga-hpo-cpu", "python").exists()


@pytest.mark.skipif(
    not _pleio_env_available(),
    reason="pga-hpo-cpu environment not available",
)
def test_extract_hpo_pleio_backend(tmp_path: Path):
    """pleio-hpo extracts the expected HPO codes from a clinical text."""
    inp = ExtractHPOInput(
        text="The patient presented with chronic kidney disease and polyuria.",
        output_dir=tmp_path / "hpo",
        backend="pleio-hpo",
    )
    result = extract_hpo(inp)

    assert result.ok
    assert result.hpo_terms_json.exists()
    # Expected codes
    assert "HP:0012622" in result.hpo_ids  # Chronic kidney disease
    assert "HP:0000103" in result.hpo_ids  # Polyuria


@pytest.mark.skipif(
    not _pleio_env_available(),
    reason="pga-hpo-cpu environment not available",
)
def test_extract_hpo_pleio_multiple_terms(tmp_path: Path):
    """Multiple HPO terms are extracted in one pass."""
    inp = ExtractHPOInput(
        text="Polyuria, polydipsia, and chronic kidney disease.",
        output_dir=tmp_path / "hpo",
        backend="pleio-hpo",
    )
    result = extract_hpo(inp)

    assert result.ok
    assert len(result.hpo_ids) >= 3
    # No duplicates
    assert len(result.hpo_ids) == len(set(result.hpo_ids))
