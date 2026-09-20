# Demo NPHP1: hallazgo del experimento lineal

Documento de resultados del primer experimento end-to-end del sistema.
Es un resultado negativo, y eso es importante: fija por qué el
consenso lineal no basta y justifica el grafo local.

## 1. Objetivo del experimento

Demostrar que la pangenomización local lineal **recupera variantes**
que un pipeline lineal contra GRCh38 pierde por reference bias.

Fixture: NPHP1, chr2:110,000,000-110,025,000, con:
- Cohorte simulada con frecuencias reales de gnomAD MID.
- Pedigree de 3 generaciones, 10 individuos.
- Individuo C2 con la variante 2-110008979-A-T forzada (AF_MID=4.4%, 0/1).
- BAM sintético de C2 con cobertura uniforme 30x.

## 2. Pipeline ejecutado

    fetch_roi          -> BAM de C2 ya disponible en el fixture
    local_pangenome    -> consenso enriquecido del ROI (24 SNVs MENA aplicados)
    call_variants      -> contra GRCh38
    call_variants      -> contra enriched_roi.fa
    compare_vcfs       -> delta categorizado

Referencias:
- GRCh38: 25000 bases.
- Enriched: 25000 bases (tras excluir indels que desplazarían coordenadas).
- 24 SNVs aplicados, 0 desplazamiento.

## 3. Resultado

Calling contra GRCh38: 18 variantes.

Calling contra enriched: 15 variantes.

compare_vcfs:
    recovered:   14
    lost:        17
    consistent:   1
    ------------------
    total:       32

## 4. Análisis de las categorías

### 4.1 La variante objetivo (8979) es `consistent`

La variante 2-110008979-A-T aparece en **ambos callings** con el
mismo genotipo 0/1. No hay ganancia: el caller contra GRCh38 la
detecta sin dificultad. La variante no es "difícil".

### 4.2 Las 14 "recovered" son falsos positivos

Caso por caso, todas las variantes marcadas como "recovered"
corresponden a posiciones donde:

- El consenso enriquecido tiene un alelo MENA común.
- C2 tiene el alelo GRCh38 en su BAM.
- El caller contra enriched interpreta la diferencia como variante.

Ejemplo (posición 10017):
    GRCh38:   G
    Enriched: A
    Cohort:   C2 = 0/0 en todas las muestras
    Reads:    todas match contra GRCh38
    Calling contra enriched reporta: "recovered G→A, GT=1/1"

C2 no porta la variante. El "recovered" es un artefacto.

### 4.3 Las 17 "lost" son en su mayoría artefactos simétricos

Misma posición, alelos invertidos. Ejemplo (posición 7032):
    GRCh38:   C
    Enriched: T
    Reads:    mayoría T (C2 es 0/1 en el cohort)
    Calling contra GRCh38 reporta: "C→T, 0/1"
    Calling contra enriched reporta: "T→C, 0/1"

Ambas detecciones describen la misma realidad (C2 heterocigoto C/T),
pero el comparador las trata como variantes distintas porque
(REF, ALT) están invertidos.

## 5. Diagnóstico

El consenso enriquecido lineal **no es un haplotipo real**. Es una
mezcla forzada de alelos comunes de la cohorte. Cuando C2 no coincide
con esa mezcla en alguna posición, el caller contra enriched genera
falsos positivos sistemáticos.

Es exactamente el problema del "Frankenstein" que ya estaba
documentado en `docs/local_pangenome.md` sección 5, pero ahora
**demostrado empíricamente**.

## 6. Consecuencias

### 6.1 El consenso lineal no sirve para diagnosticar

No es una limitación menor: es un fallo estructural del enfoque.
Llamar variantes contra un consenso lineal forzado **produce más
falsos positivos que hallazgos reales**.

### 6.2 La solución correcta es el grafo local

Un grafo pangenómico (`vg`) representa múltiples haplotipos reales
sin inventarse uno que no existe. Es la solución conocida en la
literatura. La contribución de este proyecto es demostrar que:

- **Hace falta** (este experimento lo prueba).
- **Se puede hacer on-the-fly en un portátil** porque el problema
  está confinado al ROI.

### 6.3 El consenso lineal se mantiene como modo de triaje

`local_pangenome` sigue produciendo el consenso lineal porque es
rápido y barato. Sirve como **modo de triaje**, no como referencia
diagnóstica. El modo grafo será el modo de referencia.

## 7. Trabajo futuro: grafo local

Pipeline propuesto:

    local_pangenome    -> construye subgrafo del ROI con vg construct
    align_to_graph     -> vg giraffe: alinea lecturas del paciente al grafo
    call_from_graph    -> vg call: variantes en el marco del grafo
    compare_vcfs       -> delta con verdadero valor clínico

Ventajas del grafo:
- Representa haplotipos reales, no mezclas.
- Elimina los falsos positivos del consenso lineal.
- Mantiene la eficiencia quirúrgica (todo confinado al ROI).

Requisitos:
- Entorno `pga-pangenome` con `vg`.
- Nuevas tools: `build_local_graph`, `align_to_graph`, `call_from_graph`.
- Realineamiento contra el grafo (bwa ya no aplica).

## 8. Mensaje metodológico

El experimento produce una contribución clara:

> Los pipelines lineales tienen reference bias en poblaciones
> consanguíneas. Llamar contra un consenso lineal enriquecido no
> resuelve el problema: genera falsos positivos sistemáticos. La
> solución correcta son grafos locales del ROI, y se pueden
> construir on-the-fly en hardware modesto.

Esto posiciona el proyecto como una herramienta metodológica
honesta, no como un pipeline más.

## 9. Reproducibilidad

Los comandos exactos para reproducir este experimento están en
`docs/demo_nphp1_commands.sh`. Los hashes de los ficheros usados
están en el `enrichment_report.json` correspondiente.

Resultado esperado al reproducir:
- GRCh38 calling: 18 variantes.
- Enriched calling: 15 variantes.
- compare_vcfs: recovered=14, lost=17, consistent=1.
- Análisis caso por caso: las 14 recovered son falsos positivos.
