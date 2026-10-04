# Demo: real trio with Osteopetrosis data (hg19, chr8 subset)

First run of PublicGenomicAgent on **real (non-synthetic) BAMs**,
using the public Galaxy Training trio for autosomal recessive
Osteopetrosis.

**Source**: https://zenodo.org/records/3243160

## Data

| File | Size | Sample |
|---|---|---|
| mapped_reads_father.bam | 337 MB | father (unaffected) |
| mapped_reads_mother.bam | 296 MB | mother (unaffected) |
| mapped_reads_proband.bam | 392 MB | proband (affected) |
| Pedigree.txt | 68 B | trio |

The BAMs are exome subsets aligned to **hg19**, covering only
**chr8** (146,364,022 bp). The causal gene for Osteopetrosis
(TCIRG1) is on chr11, so the variant is **not in these BAMs**.
This demo validates the pipeline on real data, not the diagnosis.

## Reference

```bash
mkdir -p data/ref && cd data/ref
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr8.fa.gz
gunzip chr8.fa.gz
samtools faidx chr8.fa
```

## Run

```bash
pga run-trio \
    --father data/osteopetrosis/mapped_reads_father.bam \
    --mother data/osteopetrosis/mapped_reads_mother.bam \
    --proband data/osteopetrosis/mapped_reads_proband.bam \
    --reference data/ref/chr8.fa \
    --region chr8:1000000-2000000 \
    --case-id OSTEO-TRIO --label OSTEO_ROI_1 \
    --clinvar ~/clinvar.vcf.gz
```

## Result

**Runtime: 19.6 seconds** on a laptop (with ClinVar annotation).

13 tool calls, all ok:

- qc_bam x3
- fetch_roi x3
- call_variants x3
- joint_call x1
- annotate_variants x1
- mendelian_filter x1
- prioritize_variants x1

The pipeline detected **3 autosomal recessive candidate
variants** with the canonical pattern:

| Locus | Change | Father | Mother | Proband | Score |
|---|---|---|---|---|---|
| chr8:1905823 | A>G | 0/1 | 0/1 | 1/1 | 3.0 |
| chr8:1905132 | G>A | 0/1 | 0/1 | 1/1 | 3.0 |
| chr8:1729582 | CTGGTGG>CTGGTGGTGG | 0/1 | 0/1 | 1/1 | 3.0 |

All three are classified as `auto_rec_hom` by `mendelian_filter`.
None are in ClinVar, so their ClinVar significance is empty.

## What this proves

- PublicGenomicAgent runs end-to-end on real BAMs.
- The pipeline scales to 350 MB per sample for a 1 Mb ROI in 20 s.
- Mendelian segregation works on real variant calls.
- ClinVar annotation is propagated from the annotated joint VCF
  to the auto_rec_hom VCF (header check).

## What this does not prove

- **Diagnosis**: the causal TCIRG1 variant is not in the data.
- **Pangenome**: no graph was built for this run.

## Next steps

- Download a real trio with full coverage of TCIRG1 (chr11, hg19).
- Benchmark linear vs pangenome on the same trio.
- Add `--clinvar` to the report to show any Pathogenic hits.

