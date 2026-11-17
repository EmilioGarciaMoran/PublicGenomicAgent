"""Benchmark justo lineal vs pangenoma local: la variante objetivo nunca está en el grafo."""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest

from publicgenomicagent.env.micromamba import bin_path

pytestmark = pytest.mark.integration
ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "scripts"))


def test_divergent_haplotype_recovered_only_with_local_graph(tmp_path):
    if not (bin_path("pga-pangenome", "vg").exists() and bin_path("pga-hts", "wgsim").exists()):
        pytest.skip("vg o wgsim no disponibles")
    import benchmark_local_pangenome as bm
    out = tmp_path / "r.json"
    assert bm.main(["--workdir", str(tmp_path / "w"), "--divergent-spacing", "12",
                    "--seed", "1", "--json", str(out)]) == 0
    r = json.loads(out.read_text())
    # la respuesta no está en el grafo: ninguna posición objetivo en el VCF de la cohorte
    with gzip.open(tmp_path / "w" / "cohort.vcf.gz", "rt") as f:
        cohort_pos = {int(l.split("\t")[1]) for l in f if not l.startswith("#")}
    lin, pan = r["branches"]["linear"], r["branches"]["local_pangenome"]
    assert not {lin[k]["pos"] for k in ("control", "divergent", "duplication")} & cohort_pos
    assert lin["control"]["called"] and pan["control"]["called"]
    assert not lin["divergent"]["called"] and pan["divergent"]["called_mapq20"]
    assert pan["false_positives"] <= lin["false_positives"]
