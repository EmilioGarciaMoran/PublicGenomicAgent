#!/usr/bin/env bash
# Reproduce docs/demo_nphp1_trio.md.
set -euo pipefail

PGA="$HOME/.pga/envs/pga-core/bin/pga"
SAMTOOLS="$HOME/.pga/envs/pga-hts/bin/samtools"
BCFTOOLS="$HOME/.pga/envs/pga-hts/bin/bcftools"
VG="$HOME/.pga/envs/pga-pangenome/bin/vg"
WORK="/tmp/pga_demo_trio"
REPO="$HOME/PublicGenomicAgent"

mkdir -p "$WORK/cohort"
cd "$WORK/cohort"

echo "=== 1. Fixture ==="
~/.pga/envs/pga-core/bin/python "$REPO/tests/fixtures/make_cohort_fixture.py" "$WORK/cohort"

echo "=== 2. Referencia ==="
~/.pga/envs/pga-core/bin/python - << 'PYEOF'
from pathlib import Path
src = Path.home() / ".pga/cache/reference/hg38_chr2_110000001_110025000.fa"
dst = Path("/tmp/pga_demo_trio/cohort/nphp1_ref.fa")
lines = src.read_text().splitlines()
seq = "".join(l for l in lines if not l.startswith(">"))
with dst.open("w") as f:
    f.write(">chr2_roi\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")
PYEOF
$SAMTOOLS faidx "$WORK/cohort/nphp1_ref.fa"

echo "=== 3. BAMs del trío ==="
for s in C2 F1 M1; do
  sl=$(echo "$s" | tr '[:upper:]' '[:lower:]')
  ~/.pga/envs/pga-core/bin/python "$REPO/tests/fixtures/make_patient_bam.py" \
    "$WORK/cohort/cohort.vcf.gz" "$WORK/cohort/nphp1_ref.fa" \
    "$s" "$WORK/cohort/${sl}.bam" 30
  $SAMTOOLS fastq -@ 2 "$WORK/cohort/${sl}.bam" > "$WORK/cohort/${sl}.fastq" 2>/dev/null
done

echo "=== 4. Grafo ==="
$PGA tool build-local-graph \
  --reference "$WORK/cohort/nphp1_ref.fa" \
  --cohort "$WORK/cohort/cohort.vcf.gz" \
  --out "$WORK/roi_tool.vg"

echo "=== 5. Alinear cada individuo ==="
for s in c2 f1 m1; do
  $PGA tool align-to-graph \
    --graph "$WORK/roi_tool.vg" \
    --reads "$WORK/cohort/$s.fastq" \
    --out "$WORK/${s}_trio.gam"
done

echo "=== 6. Pack + call por individuo ==="
for s in c2 f1 m1; do
  su=$(echo "$s" | tr '[:lower:]' '[:upper:]')
  $VG pack -x "$WORK/roi_tool.xg" -g "$WORK/${s}_trio.gam" -Q 5 -o "$WORK/${s}.pack"
  $VG call "$WORK/roi_tool.xg" -k "$WORK/${s}.pack" > "$WORK/${s}_trio.vcf" 2>/dev/null
  sed "s/^\(#CHROM.*\)\tSAMPLE\$/\1\t${su}/" "$WORK/${s}_trio.vcf" > "$WORK/${s}_named.vcf"
  $BCFTOOLS view -Oz -o "$WORK/${s}_named.vcf.gz" "$WORK/${s}_named.vcf"
  $BCFTOOLS index -t "$WORK/${s}_named.vcf.gz"
done

echo "=== 7. Merge ==="
$BCFTOOLS merge -m all \
  "$WORK/c2_named.vcf.gz" "$WORK/f1_named.vcf.gz" "$WORK/m1_named.vcf.gz" \
  > "$WORK/trio_merged.vcf"

echo "=== 8. Verificar herencia ==="
echo "Variantes por individuo:"
for s in c2 f1 m1; do
  su=$(echo "$s" | tr '[:lower:]' '[:upper:]')
  echo "  $su: $(grep -vc '^#' "$WORK/${s}_trio.vcf")"
done
echo "  Tras merge: $(grep -vc '^#' "$WORK/trio_merged.vcf")"

echo
echo "Variante objetivo 8979:"
grep -v "^##" "$WORK/trio_merged.vcf" | grep -P "\t8979\t"

echo
echo "Resultado esperado:"
echo "  C2: 26  F1: 29  M1: 23  merge: 31"
echo "  8979: C2=0/1  F1=0/1  M1=./."
