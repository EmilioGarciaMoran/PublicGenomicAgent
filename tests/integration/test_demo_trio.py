"""Test de integración end-to-end sobre el trío sintético.

Este test es el blindaje del pipeline completo:
  - Genera el trío sintético con make_demo_trio.py
  - Ejecuta AgentLoop con RuleBasedPlanner sobre el estado real
  - Verifica que el pipeline produce auto_rec_hom con el
    variant recesivo correcto (0/1, 0/1, 1/1)

Requiere entornos micromamba reales. Marcado como `integration`.
En CI rápido:
    pytest -m "not integration"
En local:
    pytest tests/integration/ -v
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from publicgenomicagent.agent.loop import AgentLoop
from publicgenomicagent.agent.planner import RuleBasedPlanner
from publicgenomicagent.agent.state import AgentState, SessionContext
from publicgenomicagent.env.paths import registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.tools.base import CaseManifest


pytestmark = pytest.mark.integration


EXPECTED_TOOL_SEQUENCE = [
    "qc_bam", "qc_bam", "qc_bam",
    "fetch_roi", "fetch_roi", "fetch_roi",
    "call_variants", "call_variants", "call_variants",
    "joint_call",
    "mendelian_filter",
]


def _run_pipeline(demo_trio: dict, tmp_path: Path) -> tuple[AgentState, SessionContext]:
    """Construye el estado y ejecuta el loop determinista."""
    case = CaseManifest.model_validate(demo_trio["case"])
    state = AgentState(case=case)
    state.bams = {
        "father": demo_trio["father_bam"],
        "mother": demo_trio["mother_bam"],
        "proband": demo_trio["proband_bam"],
    }
    state.references["hg38"] = demo_trio["ref_fa"]

    # El ROI del loop se resuelve a partir de case.candidate_rois.
    # El directorio de sub-BAMs lo controlamos vía artifacts["roi_dir"]
    # para que todo quede dentro de tmp_path.
    roi_dir = tmp_path / "roi"
    state.artifacts["roi_dir"] = roi_dir

    # El directorio de salida de mendelian también lo forzamos.
    mendelian_dir = tmp_path / "mendelian"

    runtime = ToolRuntime(load_registry(registry_path()))
    ctx = SessionContext(state=state, runtime=runtime)

    planner = RuleBasedPlanner(roi_dir=roi_dir)
    # Sobreescribimos el _rule_mendelian para usar nuestro output_dir.
    # (o simplemente dejamos que use el default y comprobamos allí).
    loop = AgentLoop(planner=planner, max_steps=30)
    loop.run(ctx)

    return state, ctx


def _read_vcf_genotypes(vcf: Path) -> list[dict]:
    """Lee un VCF con bcftools y devuelve lista de variantes con GT."""
    import subprocess
    from publicgenomicagent.env.micromamba import bin_path

    bcftools = bin_path("pga-hts", "bcftools")
    samples = subprocess.run(
        [str(bcftools), "query", "-l", str(vcf)],
        capture_output=True, text=True, check=True,
    ).stdout.strip().splitlines()

    body = subprocess.run(
        [str(bcftools), "view", "-H", str(vcf)],
        capture_output=True, text=True, check=True,
    ).stdout.strip()

    variants = []
    for line in body.splitlines():
        f = line.split("\t")
        fmt = f[8].split(":")
        gt_idx = fmt.index("GT")
        genotypes = {}
        for i, s in enumerate(samples):
            gt = f[9 + i].split(":")[gt_idx]
            genotypes[s] = gt
        variants.append({
            "chrom": f[0],
            "pos": int(f[1]),
            "ref": f[3],
            "alt": f[4],
            "genotypes": genotypes,
        })
    return variants


def test_demo_trio_runs_full_pipeline(demo_trio: dict, tmp_path: Path) -> None:
    """El pipeline completo se ejecuta sin errores."""
    state, ctx = _run_pipeline(demo_trio, tmp_path)

    # Todas las tool_calls deben ser ok
    failed = [c for c in state.tool_calls if not c.ok]
    assert failed == [], f"tool_calls fallidas: {[c.tool_name for c in failed]}"

    # La secuencia de tools debe ser la esperada
    actual = [c.tool_name for c in state.tool_calls]
    assert actual == EXPECTED_TOOL_SEQUENCE, (
        f"secuencia inesperada:\n  esperada: {EXPECTED_TOOL_SEQUENCE}\n"
        f"  actual:   {actual}"
    )


def test_demo_trio_produces_joint_vcf(demo_trio: dict, tmp_path: Path) -> None:
    """El joint_call genera un VCF con las 3 muestras."""
    state, _ = _run_pipeline(demo_trio, tmp_path)

    assert "trio" in state.vcfs, f"vcfs keys: {list(state.vcfs.keys())}"
    trio_vcf = state.vcfs["trio"]
    assert trio_vcf.exists(), f"no existe {trio_vcf}"

    variants = _read_vcf_genotypes(trio_vcf)
    assert len(variants) == 1, f"esperaba 1 variante, hay {len(variants)}"

    v = variants[0]
    assert v["chrom"] == "chr1"
    assert v["pos"] == 1100
    assert v["ref"] == "A"
    assert v["alt"] == "G"
    assert v["genotypes"] == {
        "father": "0/1",
        "mother": "0/1",
        "proband": "1/1",
    }


def test_demo_trio_classifies_auto_rec_hom(demo_trio: dict, tmp_path: Path) -> None:
    """mendelian_filter clasifica el variant como auto_rec_hom."""
    state, _ = _run_pipeline(demo_trio, tmp_path)

    # La salida de mendelian_filter queda registrada en state.outputs
    mf_out = state.latest("mendelian_filter")
    assert mf_out is not None

    counts = mf_out.get("counts", {})
    assert counts.get("auto_rec_hom") == 1, f"counts: {counts}"
    assert counts.get("de_novo") == 0
    assert counts.get("auto_dom") == 0

    # El VCF auto_rec_hom contiene la variante
    auto_rec_hom_vcf = Path(mf_out["auto_rec_hom_vcf"])
    assert auto_rec_hom_vcf.exists()

    variants = _read_vcf_genotypes(auto_rec_hom_vcf)
    assert len(variants) == 1
    assert variants[0]["pos"] == 1100
    assert variants[0]["genotypes"]["proband"] == "1/1"
    assert variants[0]["genotypes"]["father"] == "0/1"
    assert variants[0]["genotypes"]["mother"] == "0/1"


def test_demo_trio_session_is_serializable(demo_trio: dict, tmp_path: Path) -> None:
    """El estado se puede volcar a JSON y recargar sin pérdida."""
    state, _ = _run_pipeline(demo_trio, tmp_path)

    session_path = tmp_path / "session.json"
    state.to_json(session_path)
    assert session_path.exists()

    data = json.loads(session_path.read_text())
    assert data["case"]["case_id"] == "DEMO_TRIO"
    assert len(data["tool_calls"]) == 11
    # Todas las tool_calls ok
    assert all(c["ok"] for c in data["tool_calls"])
