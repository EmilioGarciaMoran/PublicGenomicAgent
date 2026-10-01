# Quickstart

A reproducible end-to-end walkthrough of PublicGenomicAgent on
a synthetic trio with a recessive variant. No real patient data
required. Runs in about 30 seconds on a laptop.

For the full policy on privacy and network use, see
[`docs/privacy.md`](privacy.md).

---

## What you will do

1. Install the base environment (`pga-core`).
2. Bootstrap the tool environments (`pga-hts`, etc.).
3. Verify the installation.
4. Generate a synthetic trio with a recessive variant.
5. Run the full analysis pipeline with a single command.
6. Inspect the results (VCFs, session trace).

## Prerequisites

| Requirement | Notes |
|---|---|
| Linux or macOS | Tested on Ubuntu and macOS |
| `micromamba` in PATH | https://mamba.readthedocs.io/ |
| Python 3.11+ | Provided by `pga-core` |
| (Optional) Ollama | Only for `--planner llm` |

## 2. Install

From the repository root:

```bash
./setup.sh
```

This creates `~/.pga/envs/pga-core`, installs the package in
editable mode, and bootstraps `pga-hts`, `pga-pangenome` and
`pga-mendelian`. It is idempotent: re-running it is safe.

If you prefer manual steps:

```bash
micromamba create -y -p ~/.pga/envs/pga-core -f envs/pga-core.yml
micromamba run -p ~/.pga/envs/pga-core pip install -e ".[dev]"
~/.pga/envs/pga-core/bin/pga env bootstrap pga-hts
```

## 3. Verify

```bash
~/.pga/envs/pga-core/bin/pga env verify
```

You should see `pga-core`, `pga-hts`, `pga-pangenome` and
`pga-mendelian` marked as `ok`.

Optional: check the local LLM (requires Ollama running with a
model already pulled, e.g. `ollama pull qwen2.5:3b`):

```bash
~/.pga/envs/pga-core/bin/pga llm ping
```

## 4. Generate the synthetic trio

```bash
python3 scripts/make_demo_trio.py
```

This creates the following files under `results/demo_trio/`:

| File | Contents |
|---|---|
| `ref.fa` | Synthetic reference (`chr1:1-10000`, all `A`) |
| `father.bam` + `.bai` | BAM with variant at `0/1` |
| `mother.bam` + `.bai` | BAM with variant at `0/1` |
| `proband.bam` + `.bai` | BAM with variant at `1/1` |
| `case.json` | `CaseManifest` of the family |

The variant is a single SNV at `chr1:1100` (`A>G`), inside the
ROI `chr1:1000-1200`. The inheritance pattern is the canonical
autosomal recessive case: both parents are healthy carriers,
the affected proband is homozygous for the alternate allele.

## 5. Run the pipeline

A single command runs the entire analysis:

```bash
~/.pga/envs/pga-core/bin/pga run-trio \\
    --father results/demo_trio/father.bam \\
    --mother results/demo_trio/mother.bam \\
    --proband results/demo_trio/proband.bam \\
    --reference results/demo_trio/ref.fa \\
    --region chr1:1000-1200 \\
    --case-id DEMO_TRIO --label DEMO1
```

The loop executes **11 automatic steps**:

```
qc_bam              x 3   (one per BAM)
fetch_roi           x 3   (one per sample, over the ROI)
call_variants       x 3   (one per sub-BAM)
joint_call          x 1   (multi-sample VCF)
mendelian_filter    x 1   (segregation analysis)
```

The final state is saved to `results/DEMO_TRIO.session.json`.

## 6. Inspect the results

The joint VCF (multi-sample):

```bash
bcftools view -H results/roi/joint_chr1_1000_1200.vcf.gz
```

Expected:

```
chr1  1100  .  A  G  ...  GT:PL:DP:AD
  0/1:...:63:23,40    # father
  0/1:...:63:23,40    # mother
  1/1:...:63:0,63     # proband
```

The recessive candidates:

```bash
bcftools view -H results/mendelian/auto_rec_hom.vcf.gz
```

This file contains the variant `chr1:1100 A>G` classified as
autosomal recessive homozygous. The other three VCFs
(`de_novo`, `auto_dom`, `x_linked_rec`) are empty, as expected.

Inspect the full tool-call trace:

```bash
python3 -c "import json; d=json.load(open('results/DEMO_TRIO.session.json')); \\
[print(f'{i}: {c[chr(34)]}{c[chr(34)]}') for i, c in enumerate(d['tool_calls'])]"
```

## 7. Next steps

- **Run with the LLM planner** (hybrid: rules first, LLM only
  when rules are exhausted):

  ```bash
  pga run-trio ... --planner llm
  ```

- **Inspect the state of a session** (dry-run of the LLM prompt):

  ```bash
  pga plan results/DEMO_TRIO.session.json --dry-run
  ```

- **Use your own data**: replace `results/demo_trio/*.bam` with
  your own BAMs and the `case.json` with your own
  `CaseManifest` (or a `.fam` file).

