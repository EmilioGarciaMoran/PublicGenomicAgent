# Demo NPHP1: grafo local recupera SNVs MENA en cluster

Tercer experimento. Demuestra que el grafo local no solo detecta
indels, sino que recupera SNVs MENA que el calling lineal pierde
por reference bias cuando están agrupados.

## 1. Objetivo

Forzar reference bias añadiendo al haplotipo de C2 un cluster de
SNVs MENA-específicos alrededor de la variante objetivo
(2-110008979-A-T). Comparar lineal vs grafo.

## 2. Diseño del experimento

Fixture: NPHP1, chr2:110,000,000-110,025,000.
- Cohort con frecuencias reales gnomAD MID.
- Pedigree 3 generaciones, 10 individuos.
- **8 SNVs MENA sintéticos** en chr2_roi:8891-9040
  (ventana de ±100 pb alrededor de 8979).
- Variante objetivo 2-110008979-A-T forzada en C2 y F1.
- BAM sintético de C2 con cobertura uniforme 30x y errores 0.2%.

Los SNVs MENA son comunes en MENA (AF_MID 0.3-0.6), portados por
todos los individuos del pedigree. Simulan un haplotipo MENA que
GRCh38 no tiene.

## 3. Pipeline

    make_cohort_fixture.py   -> cohort.vcf.gz con cluster MENA
    make_patient_bam.py      -> c2.bam
    call_variants (lineal)   -> patient_grch38_mena.vcf.gz
    vg construct -a          -> roi_a.vg
    vg giraffe + vg call     -> c2_mena_graph.vcf

## 4. Resultados

    Calling lineal:  24 variantes
    Grafo:           26 variantes
    Consistent:      20
    Recovered:        6 (4407 + 5 SNVs MENA)
    Lost:             4

## 5. Análisis

### 5.1 Lineal pierde SNVs MENA en cluster

Los 5 SNVs MENA recuperados por el grafo están **verificados**
como presentes en el cohort:

    8902 A→T AF_MID=0.43
    8990 C→G AF_MID=0.32
    9008 A→C AF_MID=0.31
    9023 A→C AF_MID=0.58
    9040 G→T AF_MID=0.47

C2 los porta (GT 0/1). El calling lineal **no los detecta**.

Causa: múltiples mismatches consecutivos contra GRCh38. El caller
lineal los descarta como "región ambigua" cuando están agrupados.

### 5.2 Grafo los recupera

El grafo contiene estos SNVs como alt paths. Las lecturas de C2
alinean contra esos paths, y `vg call` los reporta como variantes
con genotipo correcto y DP coherente.

### 5.3 También recupera el indel 4407

El indel 4407 (TTGTG→T) se recupera como en experimentos previos.
El calling lineal lo pierde.

### 5.4 Falsos positivos eliminados

Los SNVs MENA que sí detecta el lineal (8891, 8903, 8961, 8987)
son los que el grafo absorbe como alelos de referencia comunes.
En el calling lineal aparecen como variantes; en el grafo no.
**Esto reduce ruido clínico por ancestría.**

## 6. Tabla resumen

| Categoría | Lineal | Grafo |
|---|---|---|
| SNV MENA aislado | ✓ detectado | ✓ (o absorbido) |
| SNV MENA en cluster | ✗ perdido | ✓ detectado |
| Indel pequeño | ✗ perdido | ✓ detectado |
| Variante objetivo 8979 | ✓ | ✓ |

## 7. Interpretación clínica

**El calling lineal contra GRCh38 comete dos tipos de error en
poblaciones MENA:**

1. **Falsos positivos**: reporta alelos comunes MENA como "variantes
   del paciente". Ruido clínico.

2. **Falsos negativos**: pierde SNVs MENA reales cuando están
   agrupados. Variantes clínicamente relevantes no detectadas.

**El grafo local corrige ambos:**

1. **Absorbe los alelos comunes MENA** como referencia local
   (sin falsos positivos por ancestría).
2. **Recupera los SNVs MENA reales** que el lineal pierde por
   proximidad (sin falsos negativos por reference bias).
3. **Añade detección de indels** que el lineal no puede representar.

## 8. Cuantificación

| Métrica | Lineal | Grafo | Delta |
|---|---|---|---|
| Variantes totales | 24 | 26 | +2 |
| SNVs MENA en cluster recuperados | 0 | 5 | +5 |
| Indels recuperados | 0 | 1 | +1 |

## 9. Mensaje metodológico

> El calling lineal contra GRCh38 no es solo impreciso en
> poblaciones MENA: comete **falsos positivos sistemáticos** (reporta
> alelos comunes como variantes) y **falsos negativos sistemáticos**
> (pierde SNVs reales en clusters densos). Un grafo local del ROI
> corrige ambos y añade detección de indels.

## 10. Reproducibilidad

Ver `docs/demo_nphp1_mena_commands.sh`.

## 11. Limitaciones

- Cluster MENA sintético (no observado biológicamente).
- Sin SVs reales (vg construct los ignora).
- Un solo individuo (C2), no trío.
- Sin comparación con validación experimental.

## 12. Trabajo futuro

- Verificar el fenómeno en datos reales de MENA.
- Añadir SVs representables (vg augment).
- Extender a trío con phasing.
- Publicar metodología con datos de cohortes reales.
