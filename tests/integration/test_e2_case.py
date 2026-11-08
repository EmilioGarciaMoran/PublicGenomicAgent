"""E2 de punta a punta con sustitutos locales de UCSC y gnomAD (sin red)."""
from __future__ import annotations

import importlib.util
import json
import random
from pathlib import Path

import pytest

from publicgenomicagent.env.micromamba import bin_path

pysam = pytest.importorskip("pysam")
pytestmark = pytest.mark.integration
SCRIPT = Path(__file__).parents[2] / "scripts" / "make_e2_case.py"
START, END = 42123940, 42133941


def test_e2_case_end_to_end(tmp_path, monkeypatch):
    if not bin_path("pga-hts", "bcftools").exists():
        pytest.skip("pga-hts no disponible")
    import publicgenomicagent.knowledge.gnomad as gn
    monkeypatch.setattr(gn, "AF_CACHE", tmp_path / "gnomad_cache")
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(END - START + 1))
    fa = tmp_path / "roi.fa"
    fa.write_text(">chr22\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    swap = {"A": "G", "C": "T", "G": "A", "T": "C"}
    variants = []
    for pos in range(START + 200, END - 200, 37):
        ref = seq[pos - START]
        common = rng.random() < 0.5
        ac = 3000 if common else 1
        variants.append({"variant_id": f"22-{pos}-{ref}-{swap[ref]}", "pos": pos, "ref": ref,
                         "alt": swap[ref], "genome": None,
                         "exome": {"ac": ac, "an": 20000,
                                   "populations": [{"id": "mid", "ac": ac // 20, "an": 1000}]}})
    gj = tmp_path / "gnomad.json"
    gj.write_text(json.dumps(variants))

    spec = importlib.util.spec_from_file_location("e2", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.chdir(tmp_path)
    res = mod.main(["--workdir", str(tmp_path / "e2"), "--seq-fasta", str(fa),
                    "--gnomad-json", str(gj)])
    assert res["E2_passed"], res
    assert res["seeded_in_recessive_candidates"] and res["rarity_filter"]["removed_common"] > 0
    assert res["n_ref_mismatch"] == 0
    assert (tmp_path / "e2" / "e2_report.md").read_text().count("superado") == 1
