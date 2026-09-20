#!/usr/bin/env bash
# Reproduce el experimento de demo_nphp1_graph.md.
set -euo pipefail

PGA="$HOME/.pga/envs/pga-core/bin/pga"
SAMTOOLS="$HOME/.pga/envs/pga-hts/bin/samtools"
VG="$HOME/.pga/envs/pga-pangenome/bin/vg"
WORK="/tmp/pga_demo_graph"
REPO="$HOME/PublicGenomicAgent"

# Asume que el fixture ya está en $WORK/cohort (ver demo_nphp1_commands.sh)
cd "$WORK/cohort"

echo "=== 1. Extraer FASTQ del BAM ==="
$SAMTOOLS fastq -@ 2 c2.bam > c2.fastq 2>/dev/null

echo "=== 2. Construir grafo con alt paths ==="
$VG construct -r nphp1_ref.fa -v cohort.vcf.gz -m 32 -a > roi_a.vg

echo "=== 3. Indexar ==="
$VG index -x roi_a.xg -L roi_a.vg

echo "=== 4. Alinear con giraffe ==="
$VG giraffe -x roi_a.xg -f c2.fastq -o gam > c2.gam

echo "=== 5. Pack ==="
$VG pack -x roi_a.xg -g c2.gam -Q 5 -o c2.pack

echo "=== 6. Call ==="
$VG call roi_a.xg -k c2.pack > c2_graph_calls.vcf

echo "=== 7. Comparar con calling lineal ==="
grep -v "^#" c2_graph_calls.vcf | awk '{print $2}' | LC_ALL=C sort -n > /tmp/graph_pos.txt
$HOME/.pga/envs/pga-hts/bin/bcftools query -f '%POS\n' ../patient_grch38.vcf.gz | LC_ALL=C sort -n > /tmp/linear_pos.txt

echo "Grafo:  $(wc -l < /tmp/graph_pos.txt) variantes"
echo "Lineal: $(wc -l < /tmp/linear_pos.txt) variantes"
echo
echo "Recovered por grafo:"
LC_ALL=C comm -23 /tmp/graph_pos.txt /tmp/linear_pos.txt
echo
echo "Consistent:"
LC_ALL=C comm -12 /tmp/graph_pos.txt /tmp/linear_pos.txt | wc -l

echo
echo "Resultado esperado:"
echo "  Grafo: 19 variantes"
echo "  Lineal: 18 variantes"
echo "  Recovered: 1 (indel en 4407)"
echo "  Consistent: 18"
