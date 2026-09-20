#!/usr/bin/env bash
# Reproduce docs/demo_nphp1_mena.md.
set -euo pipefail

PGA="$HOME/.pga/envs/pga-core/bin/pga"
SAMTOOLS="$HOME/.pga/envs/pga-hts/bin/samtools"
BCFTOOLS="$HOME/.pga/envs/pga-hts/bin/bcftools"
VG="$HOME/.pga/envs/pga-pangenome/bin/vg"
WORK="/tmp/pga_demo_mena"
REPO="$HOME/PublicGenomicAgent"

mkdir -p "$WORK/cohort"
cd "$WORK/cohort"

echo "=== 1. Fixture con cluster MENA ==="
~/.pga/envs/pga-core/bin/python "$REPO/tests/fixtures/make_cohort_fixture.py" "$WORK/cohort"

echo "=== 2. Referencia del ROI (contig chr2_roi) ==="
~/.pga/envs/pga-core/bin/python - << 'PYEOF'
from pathlib import Path
src = Path.home() / ".pga/cache/reference/hg38_chr2_110000001_110025000.fa"
dst = Path("/tmp/pga_demo_mena/cohort/nphp1_ref.fa")
lines = src.read_text().splitlines()
seq = "".join(l for l in lines if not l.startswith(">"))
with dst.open("w") as f:
    f.write(">chr2_roi\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")
PYEOF
$SAMTOOLS faidx "$WORK/cohort/nphp1_ref.fa"

echo "=== 3. BAM de C2 ==="
~/.pga/envs/pga-core/bin/python "$REPO/tests/fixtures/make_patient_bam.py" \
  "$WORK/cohort/cohort.vcf.gz" \
  "$WORK/cohort/nphp1_ref.fa" \
  C2 "$WORK/cohort/c2.bam" 30

echo "=== 4. Calling lineal ==="
$PGA tool call-variants \
  --bam "$WORK/cohort/c2.bam" \
  --reference "$WORK/cohort/nphp1_ref.fa" \
  --out "$WORK/patient_grch38_mena.vcf.gz" \
  --region chr2_roi:1-25000 \
  --min-qual 20 --min-dp 10

echo "=== 5. Grafo ==="
$VG construct -r "$WORK/cohort/nphp1_ref.fa" -v "$WORK/cohort/cohort.vcf.gz" -m 32 -a > "$WORK/roi_a.vg"
$VG index -x "$WORK/roi_a.xg" -L "$WORK/roi_a.vg"
$SAMTOOLS fastq -@ 2 "$WORK/cohort/c2.bam" > "$WORK/c2.fastq" 2>/dev/null
$VG giraffe -x "$WORK/roi_a.xg" -f "$WORK/c2.fastq" -o gam > "$WORK/c2.gam" 2>/dev/null
$VG pack -x "$WORK/roi_a.xg" -g "$WORK/c2.gam" -Q 5 -o "$WORK/c2.pack"
$VG call "$WORK/roi_a.xg" -k "$WORK/c2.pack" > "$WORK/c2_mena_graph.vcf" 2>/dev/null

echo "=== 6. Comparación ==="
$BCFTOOLS query -f '%POS\n' "$WORK/patient_grch38_mena.vcf.gz" | LC_ALL=C sort -n > /tmp/l_pos.txt
grep -v "^#" "$WORK/c2_mena_graph.vcf" | awk '{print $2}' | LC_ALL=C sort -n > /tmp/g_pos.txt

echo "Lineal: $(wc -l < /tmp/l_pos.txt) variantes"
echo "Grafo:  $(wc -l < /tmp/g_pos.txt) variantes"
echo
echo "Recovered por grafo:"
LC_ALL=C comm -23 /tmp/g_pos.txt /tmp/l_pos.txt
echo
echo "Resultado esperado:"
echo "  Lineal: 24  Grafo: 26"
echo "  Recovered: 4407 + 5 SNVs MENA"
