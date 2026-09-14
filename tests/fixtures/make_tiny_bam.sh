#!/usr/bin/env bash
# Genera un BAM sintético mínimo con lecturas en chr1 para tests.
# Uso: make_tiny_bam.sh <samtools> <out.bam>
set -euo pipefail

SAMTOOLS="${1:?falta ruta a samtools}"
OUT="${2:?falta ruta de salida}"

SAM=$(mktemp --suffix=.sam)
trap 'rm -f "$SAM"' EXIT

{
  printf '@HD\tVN:1.6\tSO:coordinate\n'
  printf '@SQ\tSN:chr1\tLN:10000\n'
  # 5 lecturas en 100-200
  printf 'r1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n'
  printf 'r2\t0\tchr1\t120\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n'
  printf 'r3\t0\tchr1\t150\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n'
  printf 'r4\t0\tchr1\t180\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n'
  printf 'r5\t0\tchr1\t200\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n'
  # 3 lecturas en 500-600
  printf 'r6\t0\tchr1\t500\t60\t10M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII\n'
  printf 'r7\t0\tchr1\t550\t60\t10M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII\n'
  printf 'r8\t0\tchr1\t600\t60\t10M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII\n'
  # 2 lecturas en 9000-9100
  printf 'r9\t0\tchr1\t9000\t60\t10M\t*\t0\t0\tGGGGGGGGGG\tIIIIIIIIII\n'
  printf 'r10\t0\tchr1\t9100\t60\t10M\t*\t0\t0\tGGGGGGGGGG\tIIIIIIIIII\n'
} > "$SAM"

"$SAMTOOLS" view -b -o "$OUT" "$SAM"
"$SAMTOOLS" index "$OUT"
echo "OK: $OUT"
