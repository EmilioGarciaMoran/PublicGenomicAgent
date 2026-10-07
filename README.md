# PublicGenomicAgent

**A clinical-grade, privacy-first conversational agent for genomic analysis.**

Given a clinical phenotype, a family structure, and pre-aligned BAMs,
PublicGenomicAgent analyzes only the relevant regions of interest (ROIs)
and returns a prioritized variant report — with Mendelian statistics,
local pangenomization, and a fully auditable execution trace.

Designed for modest hardware: laptops, external drives, no cluster.
Air-gapped capable.

---

## What this is

PublicGenomicAgent is **not a pipeline**. A pipeline produces files.
This system produces *decisions*.

- **Input**: a pre-aligned BAM + a clinical query (phenotype, family,
  diagnostic suspicion) + optionally free-text clinical notes.
- **Output**: a quantified decision over candidate variants, with
  Mendelian segregation, local pangenomization, and complete
  traceability.

The VCF is the intermediate language. The decision is the product.

## Approach: clinical, not genomic

| Aspect | Genomic pipeline | PublicGenomicAgent |
|---|---|---|
| Starting point | FASTQ / whole exome | BAM + phenotype |
| Scope | Whole genome | Clinical-driven ROI |
| Time per case | Hours | Seconds per ROI |
| Output | Exhaustive VCF | Prioritized delta VCF |
| User | Bioinformatician | Clinical geneticist |
| Validation | Global sensitivity | Diagnostic utility |

## We do not align

The system starts with **BAMs already aligned to hg38**, produced by
your existing pipeline (DRAGEN, BWA+GATK, minimap2, etc.). We assume
the heavy alignment work is already done, on a server or in the cloud.

The value is in the clinical layer: phenotype-driven interpretation
on existing BAMs, with local pangenomization to correct for population
reference bias.

## Privacy by architecture

Three guarantees, each verified by tests:

1. **No patient data leaves the machine.** The default LLM is Ollama
   running on `localhost`. Analysis tools (samtools, bcftools, vg,
   plink) run against local binaries in isolated micromamba environments.
2. **No absolute paths are sent to the LLM.** The planner prompt
   contains basenames only (`proband.bam`), never
   `/home/hospital/patient_1234.exome.bam`.
3. **No secrets in config files.** `~/.pga/config.yaml` never holds
   API keys. If a remote provider is added later, credentials go
   through environment variables.

See [`docs/privacy.md`](docs/privacy.md) for the full policy.

## Pipeline

```
Clinical text
    ↓
extract_hpo (RDMA, GPU optional)  →  HPO terms
    ↓
phenotype_ranking (LIRICAL, CPU)  →  ranked diseases
    ↓
fetch_roi → qc_bam → call_variants → joint_call → mendelian_filter
    ↓
local_pangenome → build_local_graph → align_to_graph → call_from_graph
    ↓
compare_vcfs → prioritized delta
```

Each step is a tool of the CLI. All are deterministic, auditable, and
composable. A local LLM (Ollama or llama.cpp) can orchestrate them
through `LLMPlanner`, with a deterministic `RuleBasedPlanner` fallback
when the LLM hallucinates a tool or an argument.

## Quickstart

### Install

Requires Linux or macOS, `micromamba` in PATH.

```bash
git clone https://github.com/EmilioGarciaMoran/PublicGenomicAgent
cd PublicGenomicAgent

# 1. Create the base environment (once)
micromamba create -y -p ~/.pga/envs/pga-core -f envs/pga-core.yml
micromamba run -p ~/.pga/envs/pga-core pip install -e .

# 2. Create the tool environments
~/.pga/envs/pga-core/bin/pga env bootstrap pga-hts
~/.pga/envs/pga-core/bin/pga env bootstrap pga-pangenome
~/.pga/envs/pga-core/bin/pga env bootstrap pga-mendelian

# 3. Verify
~/.pga/envs/pga-core/bin/pga env verify
```

Optional phenotype layer (requires GPU ≥ 20 GB or a dedicated server):

```bash
~/.pga/envs/pga-core/bin/pga env bootstrap pga-phenotype
~/.pga/envs/pga-core/bin/pga phenotype bootstrap
```

### Run a trio analysis end-to-end

```bash
# 1. Generate a synthetic trio with a recessive variant (chr1:1100 A>G)
python3 scripts/make_demo_trio.py

# 2. Run the full pipeline (11 automatic steps)
pga run-trio \
    --father results/demo_trio/father.bam \
    --mother results/demo_trio/mother.bam \
    --proband results/demo_trio/proband.bam \
    --reference results/demo_trio/ref.fa \
    --region chr1:1000-1200 \
    --case-id DEMO_TRIO --label DEMO1
```

Output: `results/DEMO_TRIO.session.json` with 11 traced tool calls,
including `joint_call` and `mendelian_filter`. The variant is classified
as `auto_rec_hom` (autosomal recessive homozygous) with genotypes
`father=0/1`, `mother=0/1`, `proband=1/1`.

For a longer walkthrough, see [`docs/quickstart.md`](docs/quickstart.md).

## Tools

13 tools registered in `TOOL_REGISTRY`, all with Pydantic input/output
and deterministic execution:

| Tool | Description | Runtime |
|---|---|---|
| `fetch_roi` | Slice an indexed BAM by genomic region | yes |
| `qc_bam` | Structural QC: header, index, sort order, samples, counts | yes |
| `call_variants` | Single-sample variant calling (bcftools mpileup+call) | yes |
| `joint_call` | Multi-sample variant calling for families | yes |
| `mendelian_filter` | De novo / recessive homozygous / dominant / X-linked | no |
| `plink_validate` | PLINK-based family QC (Mendel errors, IBD) | no |
| `compare_vcfs` | Diff two VCFs with sensitivity/precision metrics | no |
| `local_pangenome` | Enrich a local reference with cohort alleles | no |
| `build_local_graph` | Build a vg graph from reference + cohort VCF | yes |
| `align_to_graph` | Align reads against a vg graph (vg giraffe) | yes |
| `call_from_graph` | Call variants from graph + pack (vg call) | yes |
| `extract_hpo` | Extract HPO terms from clinical text (RDMA) | no |
| `phenotype_ranking` | LIRICAL-based candidate ranking from HPO ± VCF | no |

List them from the CLI:

```bash
pga tool list
pga tool describe fetch_roi
```

## Virtualization

One micromamba environment per tool family. No shared binaries, no
global installs, no zombie dependencies. Manifest hashes are verified
on every bootstrap; if a manifest changes, the environment is rebuilt
from scratch.

```bash
pga env list       # show declared environments
pga env bootstrap  # create one (idempotent)
pga env verify     # check binaries and versions
```

See [`docs/virtualization.md`](docs/virtualization.md).

## Tests

```bash
pytest tests/unit/         # 65 unit tests, no environment needed
pytest tests/integration/  # 30 integration tests, requires micromamba
```

The integration test `tests/integration/test_demo_trio.py` runs the
full pipeline on a synthetic recessive trio and verifies the exact
`auto_rec_hom` classification. It is the regression net for the whole
system.

## Project structure

```
envs/         manifests and registry of micromamba environments
src/          publicgenomicagent Python package
scripts/      bootstrap and demo utilities
docs/         design notes and reproducible demos
tests/        unit and integration test suites
results/      generated outputs (git-ignored)
data/         user data (git-ignored)
```

## Sandbox integration

PublicGenomicAgent consumes synthetic case catalogs produced by
Sandbox, a separate project that generates trio BAMs + manifests
from SPDI notation. The contract is explicit: PGA reads
`manifest.yaml` + BAMs + reference; the generator is opaque.

```bash
pga sandbox list ~/sandbox
pga sandbox describe ~/sandbox/genes/NPHP1/NC_000002.12_110800000_290000_
pga sandbox run ~/sandbox/genes/CYP2D6/NC_000022.11_42128940_C_T \
    --reference ~/sandbox/test_ref.fa
```

Example (CYP2D6 trio, pharmacogenomic variant):

```
CYP2D6 (NC_000022.11:42128940:C:T)
Loop: 12 iteraciones
  OK qc_bam  (x3)
  OK fetch_roi  (x3)
  OK call_variants  (x3)
  OK joint_call
  OK mendelian_filter
  OK prioritize_variants
```

Result: one recessive candidate at `chr2:42128941 A>T` with the
expected genotypes (`father 0/1`, `mother 0/1`, `proband 1/1`).

The full pipeline is covered by
`tests/integration/test_sandbox_run.py`.

## Documentation

Design notes:

- [`docs/virtualization.md`](docs/virtualization.md) — isolated environments.
- [`docs/local_pangenome.md`](docs/local_pangenome.md) — local pangenomization.
- [`docs/decision_not_pipeline.md`](docs/decision_not_pipeline.md) — why this is different.
- [`docs/mendelian_segregation.md`](docs/mendelian_segregation.md) — clinical genetics as logic.
- [`docs/phenotype_layer.md`](docs/phenotype_layer.md) — HPO extraction and ranking.
- [`docs/audit_and_training.md`](docs/audit_and_training.md) — session logging.

Usage guides:

- [`docs/installation.md`](docs/installation.md) — detailed setup and troubleshooting.
- [`docs/quickstart.md`](docs/quickstart.md) — end-to-end walkthrough.
- [`docs/use_cases.md`](docs/use_cases.md) — four clinical scenarios.
- [`docs/pitch.md`](docs/pitch.md) — **formal project pitch** (3 pillars, benchmarks, roadmap).
- [`docs/privacy.md`](docs/privacy.md) — privacy policy.
- [`docs/known_issues.md`](docs/known_issues.md) — known issues (BAQ, etc.).

Reproducible demos:

- [`docs/demo_trio_recessive.md`](docs/demo_trio_recessive.md) — synthetic recessive trio (reproducible).
- [`docs/demo_nphp1.md`](docs/demo_nphp1.md) — linear analysis on NPHP1.
- [`docs/demo_nphp1_graph.md`](docs/demo_nphp1_graph.md) — linear vs. graph.
- [`docs/demo_nphp1_mena.md`](docs/demo_nphp1_mena.md) — graph recovers MENA SNVs.
- [`docs/demo_nphp1_trio.md`](docs/demo_nphp1_trio.md) — joint calling of a trio.

## Status

**Active development. 13 tools, 5 isolated environments, 95 tests passing.**

- Genomic layer (no GPU): complete and tested end-to-end.
- LLM orchestration (Ollama local): functional, with deterministic fallback.
- Phenotype layer (HPO extraction with RDMA): requires GPU ≥ 20 GB.
  See [`docs/phenotype_layer.md`](docs/phenotype_layer.md) §14.

## License

Apache 2.0. See [`LICENSE`](LICENSE).

