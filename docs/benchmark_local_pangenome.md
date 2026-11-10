# Linear alignment vs local pangenome: a fair benchmark

Run on 2026-11-10 with `scripts/benchmark_local_pangenome.py --sweep`,
vg 1.76.1, samtools 1.19.2, bcftools 1.19. About 2.5 minutes on a laptop.
This replaces [`benchmark_duplications.md`](benchmark_duplications.md),
whose result could not be reproduced from the repository and was
probably circular (see "Why the earlier benchmark was replaced").

## Design

Both branches use the same reads, the same aligner and the same caller.
Only the graph used for alignment differs:

| | Linear | Local pangenome |
|---|---|---|
| Graph | reference only | reference + cohort alleles (VCF) |
| Aligner | `vg giraffe`, reads projected (surjected) to the linear reference as BAM | same |
| Caller | `bcftools mpileup` + `call -mv` on the BAM | same |

The target variants are **never** in the cohort VCF: the graph cannot
contain the answer. Variants are called from the reads, not genotyped
from the graph.

Synthetic 30 kb reference, homozygous sample, paired reads 2×150 bp,
30×, error 0.2 % (`wgsim`, no extra mutations). Three scenarios per run:

1. **Control**: novel SNV in a unique region.
2. **Divergent haplotype** (reference bias): novel SNV inside a 1.5 kb
   haplotype that differs from the reference by one cohort SNV every
   100, 50, 25 or 12 bp (about 1, 2, 4 and 8 %), plus two 6 bp deletions
   and two 5 bp insertions.
3. **Segmental duplication**: two 3 kb copies differing at 30 sites
   (1 %). The novel SNV is in copy 1, inside a gene-conversion tract
   where the sample carries copy-2 bases (these converted alleles are
   in the cohort VCF).

10 seeds per divergence level; each seed changes the reference, the
alleles and the reads.

## Results

| Scenario | Divergent SNV every | Branch | Called | Called (MAPQ>=20) | Mean depth | Mean MAPQ |
|---|---|---|---|---|---|---|
| control | 100 bp | linear | 10/10 | 10/10 | 29 | 60 |
| control | 100 bp | local_pangenome | 10/10 | 10/10 | 29 | 60 |
| control | 50 bp | linear | 10/10 | 10/10 | 29 | 60 |
| control | 50 bp | local_pangenome | 10/10 | 10/10 | 29 | 60 |
| control | 25 bp | linear | 10/10 | 10/10 | 29 | 60 |
| control | 25 bp | local_pangenome | 10/10 | 10/10 | 29 | 60 |
| control | 12 bp | linear | 10/10 | 10/10 | 29 | 60 |
| control | 12 bp | local_pangenome | 10/10 | 10/10 | 29 | 60 |
| divergent | 100 bp | linear | 10/10 | 10/10 | 30.6 | 59.9 |
| divergent | 100 bp | local_pangenome | 10/10 | 10/10 | 30.6 | 60 |
| divergent | 50 bp | linear | 10/10 | 10/10 | 30.6 | 56.1 |
| divergent | 50 bp | local_pangenome | 10/10 | 10/10 | 30.6 | 60 |
| divergent | 25 bp | linear | 6/10 | 1/10 | 8.7 | 8.7 |
| divergent | 25 bp | local_pangenome | 10/10 | 10/10 | 30.6 | 60 |
| divergent | 12 bp | linear | 0/10 | 0/10 | 0 | 0.0 |
| divergent | 12 bp | local_pangenome | 10/10 | 10/10 | 30.6 | 60 |
| duplication | 100 bp | linear | 1/10 | 0/10 | 1.5 | 0.4 |
| duplication | 100 bp | local_pangenome | 1/10 | 0/10 | 28.8 | 2.0 |
| duplication | 50 bp | linear | 0/10 | 0/10 | 1.2 | 0.3 |
| duplication | 50 bp | local_pangenome | 1/10 | 0/10 | 30.1 | 2.2 |
| duplication | 25 bp | linear | 1/10 | 0/10 | 1.6 | 0.5 |
| duplication | 25 bp | local_pangenome | 1/10 | 0/10 | 29 | 2.0 |
| duplication | 12 bp | linear | 0/10 | 0/10 | 1.5 | 0.3 |
| duplication | 12 bp | local_pangenome | 1/10 | 0/10 | 30.1 | 2.1 |

False positives (sum over seeds): 100 bp: linear 51, pangenome 12; 50 bp: linear 41, pangenome 8; 25 bp: linear 29, pangenome 8; 12 bp: linear 16, pangenome 1

"Called (MAPQ>=20)": called when reads with mapping quality below 20
are ignored, as many clinical pipelines do.

## Interpretation

- **Divergent haplotype.** At about 1–2 % divergence both branches call
  the variant every time. From about 4 % the linear alignment loses
  reads and mapping quality (6/10 called, 1/10 with MAPQ ≥ 20, mean
  depth 8.7) and at about 8 % it loses the region entirely (0/10). The
  local pangenome called it 10/10 at full depth and MAPQ 60 at every
  level. Typical human divergence is ~0.1 %; levels of 4–8 % occur in
  hypervariable regions such as HLA, so this is where a local pangenome
  is expected to matter.
- **Segmental duplication.** Neither branch calls the variant reliably
  (≤ 1/10, never with MAPQ ≥ 20). The pangenome places the reads in the
  right copy (depth ~30 vs ~1.5), but with mapping quality ~2, and
  bcftools caps base quality by mapping quality. A local pangenome does
  not solve segmental duplications with short reads on its own.
- **False positives.** Fewer with the pangenome at every divergence
  level (1–12 vs 16–51, summed over 10 seeds), mostly paralogous
  mismapping in the duplication.

## Limits

- Synthetic reference and reads; homozygous sample; one region type per
  scenario. Not yet tested on real data or on structural variants.
- The cohort VCF contains exactly the alleles the sample carries. With a
  real population VCF the sample may carry alleles absent from it, and
  the gain will be smaller.
- One aligner (vg giraffe) in both branches; a linear aligner such as
  bwa-mem may behave differently from giraffe on a reference-only graph.

## Why the earlier benchmark was replaced

`benchmark_duplications.md` reported that the linear pipeline found 0
variants and the pangenome 1. Its generator and demo module were never
versioned. It called variants with `vg call`, which genotypes the sites
present in the graph; without graph augmentation it can only report a
variant that is already in the graph. The reads had MAPQ 0 set by hand
rather than by an aligner.

## Reproduce

```bash
python scripts/benchmark_local_pangenome.py --workdir /tmp/pgbench --sweep
python scripts/benchmark_local_pangenome.py --workdir /tmp/pgbench1 --divergent-spacing 25 --seed 1
```
