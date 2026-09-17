# Pangenomización local del ROI

Este documento describe el diseño de `local_pangenome`, la pieza central
de PublicGenomicAgent. Todo lo demás (fetch_roi, qc, call_variants,
mendelian) existe para alimentar o explotar esta capacidad.

---

## 1. Motivación

El genoma de referencia GRCh38 está sesgado hacia poblaciones europeas.
En poblaciones del Middle East, con alta consanguinidad y alelos
regionales ausentes en GRCh38, este sesgo produce:

- Lecturas que no mapean bien a la referencia.
- Caídas de cobertura artificiales.
- Falsos negativos en el variant calling.
- Variantes reales descartadas como "ruido" o "no mapeadas".

El reference bias es especialmente severo en genes con alta diversidad
regional como NPHP1, donde GRCh38 no representa bien los haplotipos
del Middle East.

**Solución:** construir in situ un **consenso enriquecido del ROI** que
incorpore los alelos comunes de la población regional, y llamar variantes
contra ese consenso en lugar de contra GRCh38.

---

## 2. Filosofía

- **Quirúrgico**: no se construye un pangenoma global. Solo el ROI.
- **Local-first**: cabe en un portátil. Nada de grafos de TB.
- **Determinista**: dado un cohort VCF + un umbral, el consenso es único.
- **Semilla reutilizable**: el FASTA enriquecido es un artefacto primario
  que se puede re-enriquecer con otra cohorte o con otro umbral.
- **VCF como contrato**: toda métrica de valor se materializa en VCF.
  Lo que no se puede expresar en VCF no es un resultado clínico.

---

## 3. Entradas

- **Referencia base**: `hg38.fa` (o el FASTA del ROI descargado de UCSC).
- **Cohort VCF**: frecuencias poblacionales del Middle East para el ROI.
  - Puede venir de gnomAD MID, Al Mena, EGP, o cualquier cohort anotado.
  - Debe tener `INFO/AF` (o `INFO/AF_MID`) con frecuencia poblacional.
  - Opcionalmente contiene genotipos de individuos (para otros usos).
- **Región del ROI**: `chr:start-end`.
- **Umbral de AF**: por defecto `0.01` (1%). Parametrizable con `--min-af`.

---

## 4. Salidas

- **`enriched_roi.fa`** + `.fai`: la referencia enriquecida.
  - Semilla reutilizable.
  - Misma longitud que la referencia base salvo por los indels aplicados.
- **`enriched_diff.vcf.gz`** + `.tbi`: VCF con las variantes incorporadas.
  - Trazabilidad: qué alelos de la cohorte están en el consenso.
  - Incluye `AF` como INFO para auditoría.
- **`enrichment_report.json`**: metadatos de la operación.
  - Nº de variantes en el cohort VCF del ROI.
  - Nº de variantes filtradas por AF.
  - Nº de variantes aplicadas.
  - Hash del cohort VCF, hash de la referencia base, hash del resultado.

---

## 5. Algoritmo

Paso 1. Validar referencias: hg38.fa + .fai, cohort.vcf.gz + .tbi.

Paso 2. Extraer cohort VCF al ROI:

    bcftools view -r <region> cohort.vcf.gz -Oz -o cohort_roi.vcf.gz
    bcftools index -t cohort_roi.vcf.gz

Paso 3. Filtrar por AF poblacional:

    bcftools view -i 'INFO/AF >= <min_af>' cohort_roi.vcf.gz -Oz -o cohort_roi_common.vcf.gz

Paso 4. Excluir SVs simbólicos:

    bcftools view -v snps,indels cohort_roi_common.vcf.gz -Oz -o cohort_roi_common_snv.vcf.gz

Paso 5. Generar consenso enriquecido forzando todos los ALT como hom-alt:

    bcftools +setGT cohort_roi_common_snv.vcf.gz -- -t a -n 1/1 -Oz -o cohort_roi_forced.vcf.gz
    bcftools consensus -f hg38_roi.fa cohort_roi_forced.vcf.gz > enriched_roi.fa
    samtools faidx enriched_roi.fa

Paso 6. Emitir diff VCF para trazabilidad: cohort_roi_common_snv.vcf.gz
con una columna extra que documenta la decisión (kept / filtered_by_af).

Paso 7. Calcular hashes y escribir el reporte JSON.

Nota sobre el paso 5: forzar todos los genotipos a 1/1 produce un
"Frankenstein" que ningún individuo real tiene. Es precisamente lo
que buscamos: una referencia que contenga todos los alelos comunes
de la población, para que las lecturas del paciente (que probablemente
porta un subconjunto de esos alelos) mapeen mejor contra ella.

Nota sobre SVs: los SVs simbólicos (DEL, INS, DUP, INV) rompen
bcftools consensus por su REF=N y por afectar coordenadas. Se tratan
aparte en el futuro sv_detector. El consenso enriquecido maneja
solo SNVs e indels.

---

## 6. Modo grafo (fase posterior)

Cuando el caso no se resuelve con el consenso lineal, o el usuario
lo pide explícitamente, se puede construir un grafo local con `vg`:

    vg construct -r hg38_roi.fa -v cohort_roi.vcf.gz > roi_graph.vg
    vg index roi_graph.vg

El agente decide cuándo usar grafo en lugar de lineal:

- Casos sin diagnóstico tras el pipeline lineal.
- Sospecha de SV complejo.
- Petición explícita del usuario.

Modo lineal por defecto. Modo grafo bajo demanda.

---

## 7. Ganancia diagnóstica: VCF vs VCF

La métrica de oro del sistema no es el consenso. Es el VCF diferencial.

El consenso enriquecido es un medio. Lo que demuestra valor clínico es
que, llamando variantes contra el consenso enriquecido en lugar de contra
GRCh38, aparecen variantes en el VCF del paciente que el pipeline
tradicional no detecta.

### 7.1 Comparación canónica

Dado un paciente con BAM del ROI:

    patient_grch38.vcf.gz    <- call_variants contra GRCh38
    patient_pga.vcf.gz       <- call_variants contra enriched_roi.fa

Se calcula:

| Categoría    | Definición                                | Interpretación clínica                     |
|--------------|-------------------------------------------|--------------------------------------------|
| Recuperadas  | En patient_pga, no en patient_grch38      | Candidatas a falsos negativos por ref bias |
| Perdidas     | En patient_grch38, no en patient_pga      | Candidatas a falsos positivos por ref bias |
| Consistentes | En ambos, mismo genotipo                  | Sin cambio                                 |
| Discordantes | En ambos, distinto genotipo               | Alelos con mapeo ambiguo a GRCh38          |

### 7.2 Con ground truth (fixture)

Cuando existe verdad de referencia (como en nuestro fixture), se añaden:

- Sensibilidad = TP / (TP + FN)
- Precisión   = TP / (TP + FP)
- Recall por genotipo: hom-alt / het / hom-ref correctos.

Comparación directa:

| Pipeline             | Sensibilidad | Precisión | Genotipos correctos |
|----------------------|--------------|-----------|---------------------|
| GRCh38 (baseline)    | X%           | Y%        | Z%                  |
| Consenso enriquecido | X'%          | Y'%       | Z'%                 |
| Delta                | +dX          | +dY       | +dZ                 |

Este cuadro es el entregable principal del proyecto.

### 7.3 Producto: compare_vcfs

El comparador emite:

- delta.vcf.gz: VCF anotado donde cada variante lleva en INFO
  PGA_STATUS = recovered | lost | consistent | discordant.
- delta_report.json: contadores, sensibilidad, precisión, genotipos.
- delta.tsv: tabla plana para inspección rápida o publicación.

### 7.4 Materialización en el fixture

El fixture contiene una variante común en Middle East con AF_MID=4.4%
(NPHP1, posición 110008979). En el pipeline completo:

1. Se generan lecturas sintéticas del paciente C2 (a partir de su VCF).
2. Se llama contra GRCh38 -> se espera que esa variante se pierda o
   tenga genotipo erróneo.
3. Se llama contra enriched_roi.fa -> se espera que se recupere.
4. compare_vcfs cuantifica la diferencia.

El fixture debe demostrar la ganancia de forma numérica y reproducible.

### 7.5 Por qué VCF y no otra cosa

- El VCF es el estándar de facto en genómica clínica.
- Un genetista puede abrir delta.vcf.gz en IGV, VarSeq, etc.
- bcftools isec, bcftools stats, vcf-compare operan sobre VCF.
- Cualquier evaluación futura (ACMG, clasificación clínica) parte de VCF.
- El VCF es interoperable: se puede comparar con otros pipelines,
  con otros laboratorios, con otros métodos.

Cualquier métrica que no se materialice en VCF (mapabilidad, cobertura
media, entropía del grafo) es un medio, no un resultado.

---

## 8. Integración con el resto del pipeline

    fetch_roi         -> sub-BAM del paciente (ROI)
    qc_bam            -> validación
    local_pangenome   -> referencia enriquecida del ROI   <- NÚCLEO
    call_variants     -> llamar contra enriched_roi.fa
    call_variants     -> llamar contra GRCh38 (baseline, para comparar)
    compare_vcfs      -> delta.vcf.gz + métricas           <- BENCHMARK
    mendelian         -> filtros según pedigree
    annotate          -> anotación funcional
    report            -> informe HTML

local_pangenome se invoca:

- Después de identificar el ROI clínico.
- Antes de llamar variantes del paciente.
- Opcionalmente: si el primer calling no resuelve, se puede probar
  con otro umbral o con grafo.

---

## 9. Fixture

El fixture se genera con tests/fixtures/make_cohort_fixture.py:

- ROI: NPHP1, chr2:110,000,000-110,025,000.
- Cohort VCF: frecuencias reales de gnomAD MID (124 variantes).
- Genotipos: sintéticos, derivados por HWE y segregación mendeliana.
- Pedigree: 3 generaciones, 10 individuos, 2 afectados (C2, C3).
- SVs: 7 plausibles por tipo/tamaño (no entran en el consenso).
- Coordenadas: locales al ROI (contig chr2_roi, posiciones 1..25000).
- Referencia: descargada de UCSC, uppercaseada, .fai con samtools faidx.

Caveat importante: la cohorte simulada tiene solo 10 individuos, lo
cual da una resolución de AF de 5%. Por debajo, la distribución es
artificialmente binaria (0 o >=5%). Por eso el filtro se basa en
INFO/AF (poblacional), no en AF observada en la cohorte.

---

## 10. Decisiones de diseño abiertas

- Umbral por defecto: 0.01. Parametrizable. Se puede reevaluar
  cuando tengamos distribución real de un cohort grande.
- Comportamiento ante indels solapados: bcftools consensus emite
  warning y los salta. Documentar como limitación conocida.
- Modo grafo: fase posterior, bajo demanda.
- Re-enriquecimiento iterativo: el FASTA enriquecido puede ser input
  de una segunda pasada con otra cohorte. No implementado aún.
- Cohortes múltiples: fusionar gnomAD MID + Al Mena + EGP antes
  del filtrado. Pendiente.

---

## 11. Análisis incremental iterativo con diálogo

El sistema no es un pipeline de un solo disparo. Es una herramienta
que un genetista experto usa en ciclos cortos, refinando el análisis
a medida que aprende del caso. Cada iteración debe ser rápida
(segundos) y trazable (qué cambió respecto a la anterior).

### 11.1 Ciclo típico

    Genetista: "Tengo un trío con fallo renal y enanismo.
                Analiza NPHP1 en C2, F1 y M1."

    Agente:    fetch_roi -> qc_bam -> local_pangenome (umbral 0.01)
               -> call_variants x2 (GRCh38 y enriquecido) -> compare_vcfs
               Reporta: "3 recuperadas, 1 discordante, 0 perdidas."

    Genetista: "Sube el umbral a 0.05, sospecho de alelo raro."

    Agente:    Repite solo local_pangenome + callings + compare_vcfs.
               Reporta: "Ahora 1 recuperada, 1 discordante."

    Genetista: "Muéstrame el alineamiento de la discordante en IGV."

    Agente:    Genera igv-reports sobre la posición concreta.

    Genetista: "Descarta esa variante, es artefacto de homopolímero."

    Agente:    Anota como filtrada; excluye de próximos refinamientos.

    Genetista: "Amplía el ROI a NPHP1 + IFT140 + DYNC2H1."

    Agente:    Repite el ciclo sobre el ROI extendido.

### 11.2 Principios del modo iterativo

- Cada paso es idempotente y cacheable: repetir con los mismos
  parámetros no recalcula.
- El estado de la sesión se mantiene en un SessionState explícito:
  ROI actual, umbral, variantes descartadas por el usuario, hipótesis.
- El agente propone el siguiente paso, pero el genetista decide.
- Los comandos son cortos y en lenguaje natural cuando la UI del
  chat esté lista; en CLI son subcomandos tipados y auditables.
- Todas las decisiones del usuario quedan en el log de sesión.

### 11.3 Validez de una variante: criterios que el genetista valora

El sistema debe exponer, para cada variante candidata, los criterios
que un genetista usa para valorar fiabilidad:

- Cobertura en la posición (DP).
- Balance de alelos (AD, VAF).
- Calidad base (QUAL, BQ).
- Sesgo de hebra (SB, FS).
- Mapabilidad local (MQ, MQ0F).
- Concordancia entre llamadas GRCh38 y enriquecida.
- Presencia en cohortes poblacionales (AF_MID, gnomAD).
- Anotación funcional (post-annotate).
- Coherencia con el pedigree.

El agente no decide: presenta la evidencia y el genetista concluye.

### 11.4 Persistencia de sesión

Un SessionState (JSON en el directorio de resultados) contiene:

- ROI(s) analizados.
- Umbral(es) probados.
- Ficheros generados en cada iteración.
- Variantes descartadas por el usuario, con motivo.
- Variantes priorizadas por el usuario.
- Hipótesis activas.

El estado se puede guardar, reabrir, y retomar días después.

---

## 12. Validación

Tests de integración previstos:

- test_local_pangenome_produces_fasta: consenso generado, longitud
  coherente con indels aplicados.

- test_local_pangenome_diff_vcf_traces_alleles: el VCF diff contiene
  las variantes que se aplicaron, con su AF original.

- test_local_pangenome_af_threshold: variantes con AF < umbral quedan
  fuera del consenso.

- test_compare_vcfs_categorizes_correctly: clasificación
  recovered / lost / consistent / discordant sobre VCFs sintéticos.

- test_pga_recovers_known_variant (benchmarking de oro):
    1. Generar BAM sintético de C2 con la variante 2-110008979-A-T.
    2. Calling contra GRCh38 -> variante NO detectada o genotipo erróneo.
    3. Calling contra enriquecido -> variante SÍ detectada, genotipo correcto.
    4. compare_vcfs -> assert recovered >= 1.
