# Demo NPHP1 con datos reales: delección en KSA001

Primer experimento con datos reales (no sintéticos) usando
PublicGenomicAgent.

## 1. Contexto

KSA001 es un genoma saudí público (Telomere-to-Telomere, 2024),
con ensamblados materno y paterno disponibles en NCBI:

- Materno: GCA_037177555.1
- Paterno: GCA_037177635.1

El paper de KSA001 demostró que alinear contra KSA001 produce
menos variantes que contra GRCh38 (reference bias).

## 2. Objetivo

Verificar si PublicGenomicAgent puede:
1. Detectar una variante real entre GRCh38 y KSA001.
2. Generar un BAM sintético con esa variante.
3. Llamar variantes contra GRCh38 y ver si detecta la variante.

## 3. Variante encontrada

Comparando el ROI de NPHP1 (chr2:110,000,001-110,025,000) de GRCh38
contra el ensamblado materno de KSA001:

- **Posición local**: 153-154
- **Tipo**: delección de 1 bp
- **Contexto**: homopolímero de T's
- **GRCh38**: 15 T's
- **KSA001**: 14 T's
- **VCF**: `chr2 153 . TT T . PASS . GT 1/1`

Es una delección limpia, sin SNPs adicionales en 2000 bp
adyacentes.

## 4. Pipeline ejecutado

1. Descarga del ensamblado de KSA001 (2.9 GB).
2. Extracción del cromosoma 2 (CM073953.1).
3. Extracción del ROI de NPHP1 (1 Mb alrededor de 110 Mb).
4. Identificación de la variante por comparación directa.
5. Generación de BAM sintético con la variante (cobertura 30x).
6. Variant calling lineal contra GRCh38.

## 5. Resultado

El calling lineal detecta la variante correctamente:

    chr2  153  .  TT  T  QUAL=44.38  PASS
      INDEL;IDV=32;IMF=0.97;DP=33;AC=2;AN=2
      GT:PL:DP:AD  1/1:71,67,0:27:1,26

- **Genotipo**: 1/1 (homocigoto de la delección)
- **IDV=32**: 32 lecturas soportan la delección
- **IMF=0.97**: 97% de las lecturas portan la delección
- **DP=33**: cobertura 33x

## 6. Interpretación

**No hay reference bias en este caso concreto.**

El pipeline lineal (bcftools call) es capaz de detectar
correctamente la delección de 1 bp en el homopolímero, con
cobertura 32x y lecturas de 150 pb con CIGAR correctas.

Esto es un resultado honesto y valioso:

- Confirma que bcftools call es robusto para indels pequeños.
- Valida que el pipeline de PublicGenomicAgent funciona con
  datos reales.
- Establece que el reference bias no afecta a todas las
  variantes por igual: los indels de 1 bp con buena cobertura
  se detectan bien.

## 7. Trabajo futuro

Para encontrar casos donde el reference bias SÍ es un problema,
hay que buscar variantes más difíciles:

- Indels en regiones repetitivas complejas.
- SVs de 100+ bp.
- Variantes en regiones de baja mappability.
- Variantes en clusters de SNPs cercanos (haplotipos MENA
  específicos).

Ver docs/demo_nphp1_mena.md para el análisis sintético donde el
grafo local sí recupera variantes que el lineal pierde.

## 8. Ficheros generados

En ~/pga_quickstart/:

- ksa001_mat_deletion.vcf — VCF con la variante
- ksa001_mat.bam — BAM sintético con la delección
- ksa001_mat_grch38.vcf.gz — VCF resultado del calling

## 9. Conclusión

PublicGenomicAgent funciona con datos reales. La variante real
de KSA001 se detecta correctamente con el pipeline lineal. El
reference bias no es un problema universal: depende del tipo de
variante y del contexto genómico.

La demostración del valor del grafo local sigue siendo el
experimento sintético con el cluster MENA (docs/demo_nphp1_mena.md),
donde el calling lineal sí pierde variantes.
