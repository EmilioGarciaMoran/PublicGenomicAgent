#!/usr/bin/env bash
# Reproduce el experimento de demo_nphp1.md.
# Requiere: pga-core bootstrapeado, pga-hts bootstrapeado,
# fixture regenerado, referencia del ROI cacheada.

set -euo pipefail

PGA="$HOME/.pga/envs/pga-core/bin/pga"
SAMTOOLS="$HOME/.pga/envs/pga-hts/bin/samtools"
BCFTOOLS="$HOME/.pga/envs/pga-hts/bin/bcftools"
WORK="/tmp/pga_demo_nphp1"
REPO="$HOME/PublicGenomicAgent"

mkdir -p "$WORK/cohort" "$WORK/pangenome"

echo "=== 1. Generar fixture ==="
~/.pga/envs/pga-core/bin/python "$REPO/tests/fixtures/make_cohort_fixture.py" "$WORK/cohort"

echo "=== 2. Generar referencia del ROI (contig chr2_roi) ==="
~/.pga/envs/pga-core/bin/python - << 'PYEOF'
from pathlib import Path
src = Path.home() / ".pga/cache/reference/hg38_chr2_110000001_110025000.fa"
dst = Path("/tmp/pga_demo_nphp1/cohort/nphp1_ref.fa")
lines = src.read_text().splitlines()
seq = "".join(l for l in lines if not l.startswith(">"))
with dst.open("w") as f:
    f.write(">chr2_roi\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")
PYEOF
$SAMTOOLS faidx "$WORK/cohort/nphp1_ref.fa"

echo "=== 3. BAM sintético de C2 ==="
~/.pga/envs/pga-core/bin/python "$REPO/tests/fixtures/make_patient_bam.py" \
  "$WORK/cohort/cohort.vcf.gz" \
  "$WORK/cohort/nphp1_ref.fa" \
  C2 \
  "$WORK/cohort/c2.bam" \
  30

echo "=== 4. Consenso enriquecido ==="
$PGA tool local-pangenome \
  --cohort "$WORK/cohort/cohort.vcf.gz" \
  --reference "$WORK/cohort/nphp1_ref.fa" \
  --region chr2_roi:1-25000 \
  --out "$WORK/pangenome" \
  --min-af 0.05 \
  --af-field AF_MID \
  --contig chr2_roi

echo "=== 5. Calling contra GRCh38 ==="
$PGA tool call-variants \
  --bam "$WORK/cohort/c2.bam" \
  --reference "$WORK/cohort/nphp1_ref.fa" \
  --out "$WORK/patient_grch38.vcf.gz" \
  --region chr2_roi:1-25000 \
  --min-qual 20 --min-dp 10

echo "=== 6. Calling contra enriched ==="
$PGA tool call-variants \
  --bam "$WORK/cohort/c2.bam" \
  --reference "$WORK/pangenome/enriched_roi.fa" \
  --out "$WORK/patient_pga.vcf.gz" \
  --region chr2_roi:1-25000 \
  --min-qual 20 --min-dp 10

echo "=== 7. Comparación ==="
$PGA tool compare-vcfs \
  --baseline "$WORK/patient_grch38.vcf.gz" \
  --candidate "$WORK/patient_pga.vcf.gz" \
  --out "$WORK/delta.vcf.gz" \
  --report "$WORK/delta.json" \
  --tsv "$WORK/delta.tsv" \
  --sample C2

echo
echo "Resultado esperado:"
echo "  GRCh38: 18 variantes"
echo "  Enriched: 15 variantes"
echo "  delta: recovered=14 (falsos positivos), lost=17, consistent=1"
echo
echo "Ficheros en $WORK"
