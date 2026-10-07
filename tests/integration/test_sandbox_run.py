"""Integration test for `pga sandbox run` on a synthetic case.

Uses the `synthetic_sandbox` fixture to build a trio with a SNV
(A>G at chr1:500) and runs the pipeline. Verifies:

  - the pipeline executes the expected sequence of tools
  - auto_rec_hom.vcf.gz contains the expected variant
  - the genotype is father 0/1, mother 0/1, proband 1/1
"""
from __future__ import annotations

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
from publicgenomicagent.tools.base import (
    CaseManifest, GenomicRange, IndividualSpec,
)


pytestmark = pytest.mark.integration


def _build_state(sandbox: dict) -> AgentState:
    case_dir = sandbox["case_dir"]
    ref_fa = sandbox["ref_fa"]
    bams_dir = case_dir / "bams"

    manifest = CaseManifest(
        case_id="SYNTH_TEST",
        source="clinical",
        pedigree=[
            IndividualSpec(sample="father", sex="M", affected=False),
            IndividualSpec(sample="mother", sex="F", affected=False),
            IndividualSpec(
                sample="proband", sex="U", affected=True,
                father="father", mother="mother",
            ),
        ],
        proband="proband",
        affected_samples=["proband"],
        candidate_rois=[GenomicRange(
            chrom="chr1", start=400, end=600, label="SYNTH",
        )],
        candidate_genes=["SYNTH"],
        consanguinity=False,
    )
    state = AgentState(case=manifest)
    state.bams = {
        "father": bams_dir / "father.bam",
        "mother": bams_dir / "mother.bam",
        "proband": bams_dir / "proband.bam",
    }
    state.references["hg38"] = ref_fa
    return state


def _read_vcf_genotypes(vcf: Path) -> list[dict]:
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
        genotypes = {s: f[9 + i].split(":")[gt_idx] for i, s in enumerate(samples)}
        variants.append({
            "chrom": f[0], "pos": int(f[1]),
            "ref": f[3], "alt": f[4],
            "genotypes": genotypes,
        })
    return variants


def test_sandbox_run_end_to_end(synthetic_sandbox, tmp_path):
    """Full pipeline on a synthetic Sandbox case."""
    state = _build_state(synthetic_sandbox)
    runtime = ToolRuntime(load_registry(registry_path()))
    ctx = SessionContext(state=state, runtime=runtime)
    planner = RuleBasedPlanner()
    loop = AgentLoop(planner=planner, max_steps=20)
    loop.run(ctx)

    # All tool calls succeeded
    failed = [c for c in state.tool_calls if not c.ok]
    assert failed == [], f"Failed tools: {[c.tool_name for c in failed]}"

    # auto_rec_hom has exactly the expected variant
    mf = state.latest("mendelian_filter")
    assert mf is not None
    assert mf["counts"]["auto_rec_hom"] == 1

    vcf = Path(mf["auto_rec_hom_vcf"])
    variants = _read_vcf_genotypes(vcf)
    assert len(variants) == 1
    v = variants[0]
    assert v["chrom"] == "chr1"
    assert v["pos"] == 501  # 1-based; 0-based chr1:500 -> VCF chr1:501
    assert v["genotypes"] == {
        "father": "0/1",
        "mother": "0/1",
        "proband": "1/1",
    }
