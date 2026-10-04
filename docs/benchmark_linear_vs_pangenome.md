# Benchmark: linear vs pangenome (Osteopetrosis trio, chr8)

First quantitative comparison between the linear pipeline
(`bcftools call`) and the local pangenome pipeline
(`vg construct` + `vg giraffe` + `vg call`) on real BAMs.

## Setup

**Data**: Osteopetrosis trio (Galaxy Training, Zenodo 3243160)
**Region**: `chr8:1000000-2000000` (hg19)
**Reads**: 52,629 from the proband BAM
**Cohort VCF**: 127 variants from the joint trio VCF

## Pipelines

**Linear**:

```
samtools view -b BAM chr8:1000000-2000000
  | bcftools mpileup -B -f chr8.fa
  | bcftools call -mv
```

**Pangenome**:

```
vg construct -r chr8.fa -v cohort.vcf.gz -a -m 32
vg index -x local.xg local.vg
vg giraffe -Z local.vg -f proband.fastq -o gam
vg pack -x local.vg -g proband.gam -o proband.pack
vg call local.vg -k proband.pack -z local.xg
```

## Results

| Category | Count |
|---|---|
| Linear-only | 26 |
| Pangenome-only | 15 |
| Shared | 76 |
| Linear total | 102 |
| Pangenome total | 91 |
| Concordance | 76/91 = 83.5% |

## Interpretation

### Linear-only variants (26)

Most are low-quality artifacts:

- `chr8:1732338 C>A` — DP=1, QUAL=5.76 (1 read, artifact)
- `chr8:1497865 G>T` — DP=7, QUAL=7.21 (low quality)
- `chr8:1575157 G>T` — DP=14, QUAL=8.25 (low quality)
- `chr8:1765220 C>A` — DP=4, QUAL=12.59 (low quality)

The linear pipeline does not filter by depth or quality
by default, so it includes these.

### Pangenome-only variants (15)

Some are real (indels in repetitive regions), some are
alternative representations of shared variants:

- `chr8:1729582 C>CTGG` — this is the recessive variant at
  chr8:1729582 detected by the linear pipeline as
  `CTGGTGG>CTGGTGGTGG`. Different representation, same event.
- `chr8:1650536 >51625>51627 A>AC` — insertion in a
  homopolymer, represented with graph node coordinates.

## Conclusion

For this trio and region, the **pangenome does not add
diagnostic value over the linear pipeline**. However:

1. The pangenome is **more conservative**: it excludes the
   ~20 low-quality artifacts that the linear pipeline keeps.
2. The pangenome handles **indels in repetitive regions**
   differently (graph node coordinates).
3. To show a real advantage, the benchmark would need to be
   run on a **region with segmental duplications** (e.g.
   NPHP1 at chr2q13), where linear alignment fails.

## Next steps

- Run the same benchmark on chr2q13 (NPHP1), where the
  linear pipeline suffers from reference bias.
- Add quality filters to the linear pipeline for a fair
  comparison.
- Compute per-variant concordance with a truth set
  (GIAB or synthetic).

