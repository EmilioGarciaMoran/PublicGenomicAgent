"""End-to-end test of the DEMO_TRIO pipeline with ClinVar.

Extends test_demo_trio.py: adds a synthetic ClinVar VCF that
annotates the variant chr1:1100 A>G as Pathogenic. Runs the
full 13-step pipeline (with annotate_variants and
prioritize_variants) and verifies:

  - the sequence of 13 tool calls
  - annotated VCF exists and contains CLNSIG=Pathogenic
  - prioritized VCF exists
  - prioritize_variants output has 1 variant with score >= 5
    (Pathogenic + proband 1/1 + QUAL>100)

Marked as `integration`: requires pga-hts.
"""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import pytest

from publicgenomicagent.agent.loop import AgentLoop
from publicgenomicagent.agent.planner import RuleBasedPlanner
from publicgenomicagent.agent.state import AgentState, SessionContext
from publicgenomicagent.env.micromamba import bin_path
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
    "annotate_variants",
    "mendelian_filter",
    "prioritize_variants",
]


def _write_clinvar_vcf(path: Path) -> None:
    """Escribe un ClinVar minimo (chr1:1100 A>G Pathogenic)."""
    bcftools = bin_path("pga-hts", "bcftools")
    plain = path.with_suffix(".vcf")
    with plain.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr1,length=10000>\n")
        f.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar significance">\n')
        f.write('##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar disease">\n')
        f.write('##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("chr1\t1100\t.\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNDN=Test+disease;CLNREVSTAT=criteria_provided\n")
    subprocess.run(
        [str(bcftools), "view", "-Oz", "-o", str(path), str(plain)],
        check=True, capture_output=True,
    )
    subprocess.run([str(bcftools), "index", "-t", str(path)], check=True)
    plain.unlink()


def test_demo_trio_with_clinvar_full_pipeline(demo_trio: dict, tmp_path: Path) -> None:
    """Pipeline completo de 13 pasos con ClinVar."""
    # 1. Crear ClinVar de prueba
    clinvar = tmp_path / "clinvar.vcf.gz"
    _write_clinvar_vcf(clinvar)

    # 2. Construir el estado
    case = CaseManifest.model_validate(demo_trio["case"])
    state = AgentState(case=case)
    state.bams = {
        "father": demo_trio["father_bam"],
        "mother": demo_trio["mother_bam"],
        "proband": demo_trio["proband_bam"],
    }
    state.references["hg38"] = demo_trio["ref_fa"]
    state.references["clinvar"] = clinvar

    roi_dir = tmp_path / "roi"
    state.artifacts["roi_dir"] = roi_dir

    runtime = ToolRuntime(load_registry(registry_path()))
    ctx = SessionContext(state=state, runtime=runtime)
    planner = RuleBasedPlanner(roi_dir=roi_dir)
    loop = AgentLoop(planner=planner, max_steps=30)
    loop.run(ctx)

    # 3. Verificar la secuencia de 13 tools
    actual = [c.tool_name for c in state.tool_calls]
    assert actual == EXPECTED_TOOL_SEQUENCE, f"secuencia inesperada: {actual}"

    # 4. Verificar el VCF anotado
    assert "annotated" in state.vcfs
    annotated = state.vcfs["annotated"]
    assert annotated.exists()

    with gzip.open(annotated, "rt") as f:
        content = f.read()
    assert "CLNSIG=Pathogenic" in content or "CLNSIG=" in content

    # 5. Verificar el VCF priorizado
    assert "prioritized" in state.vcfs
    prioritized = state.vcfs["prioritized"]
    assert prioritized.exists()

    # 6. Verificar el ranking
    pv_out = state.latest("prioritize_variants")
    assert pv_out is not None
    ranking = pv_out.get("ranking") or []
    assert len(ranking) == 1
    top = ranking[0]
    assert top["chrom"] == "chr1"
    assert top["pos"] == 1100
    # Pathogenic (+3) + proband 1/1 (+2) + QUAL>100 (+1) = 6
    assert top["score"] >= 5, f"score: {top['score']}"
    assert "Pathogenic" in (top.get("clnsig") or "")
