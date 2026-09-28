# Quickstart

Ejemplo end-to-end con el fixture NPHP1: un trío (padre, madre, hijo)
con fallo renal y enanismo.

Duración: ~5 minutos.

Requisitos: PublicGenomicAgent instalado (ver docs/installation.md).

## Qué vamos a hacer

1. Generar el fixture NPHP1 (trío sintético + cohorte MENA).
2. Extraer HPO del texto clínico (si hay GPU) o usar HPO de ejemplo.
3. Rankear enfermedades candidatas con LIRICAL.
4. Cortar los BAMs al ROI de NPHP1.
5. Joint calling del trío contra GRCh38.
6. Aplicar filtros mendelianos.
7. (Opcional) Pangenómica local: grafo vs lineal.

## Preparación

    mkdir -p ~/pga_quickstart
    cd ~/pga_quickstart

Todos los ficheros generados en este tutorial viven aquí. Se pueden
borrar cuando termines.

## Paso 1. Generar el fixture NPHP1

El fixture se construye desde cero con frecuencias reales de gnomAD
MID (Middle Eastern) para el gen NPHP1 (chr2:110,000,000-110,025,000).

    ~/.pga/envs/pga-core/bin/python \
      ~/PublicGenomicAgent/tests/fixtures/make_cohort_fixture.py \
      ./cohort

Salida esperada:

    OK: cohort/cohort.vcf.gz
    SNVs: 132  SVs: 7  samples: 10
    Contig local: chr2_roi (posiciones 1..25000)

Genera tres ficheros:

- cohort/cohort.vcf.gz — cohort VCF con 10 individuos (pedigrí de
  3 generaciones) y frecuencias reales de gnomAD MID.
- cohort/pedigree.fam — estructura familiar para PLINK.
- cohort/expected.json — variantes conocidas (ground truth).

El trío de interés: F1 (padre), M1 (madre), C2 (hijo afectado).

## Paso 2. Preparar la referencia del ROI

El fixture usa coordenadas locales (chr2_roi) para evitar conflictos
con GRCh38. Descargamos la secuencia del ROI desde UCSC y renombramos
el contig.

    ~/.pga/envs/pga-core/bin/python -c "
    from pathlib import Path
    src = Path.home() / '.pga/cache/reference/hg38_chr2_110000001_110025000.fa'
    dst = Path('cohort/nphp1_ref.fa')
    lines = src.read_text().splitlines()
    seq = ''.join(l for l in lines if not l.startswith('>'))
    with dst.open('w') as f:
        f.write('>chr2_roi\n')
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + '\n')
    "

    ~/.pga/envs/pga-hts/bin/samtools faidx cohort/nphp1_ref.fa

Salida esperada:

    cohort/nphp1_ref.fa
    cohort/nphp1_ref.fa.fai

## Paso 3. Generar los BAMs del trío

Generamos un BAM sintético para cada individuo (F1, M1, C2) con
cobertura uniforme 30x y errores de secuenciación realistas
(0.2%, transiciones sesgadas).

    for s in C2 F1 M1; do
      sl=$(echo $s | tr '[:upper:]' '[:lower:]')
      ~/.pga/envs/pga-core/bin/python \
        ~/PublicGenomicAgent/tests/fixtures/make_patient_bam.py \
        cohort/cohort.vcf.gz cohort/nphp1_ref.fa $s cohort/${sl}.bam 30
    done

Salida esperada:

    OK: /home/<user>/pga_quickstart/cohort/c2.bam
    Lecturas: 5000  Sample: C2  Contig: chr2_roi  Cobertura objetivo: 30x
    OK: /home/<user>/pga_quickstart/cohort/f1.bam
    Lecturas: 5000  Sample: F1  Contig: chr2_roi  Cobertura objetivo: 30x
    OK: /home/<user>/pga_quickstart/cohort/m1.bam
    Lecturas: 5000  Sample: M1  Contig: chr2_roi  Cobertura objetivo: 30x

## Paso 4. Extraer fenotipo del texto clínico

El genetista escribe una descripción clínica del caso. El sistema
extrae términos HPO usando RDMA (requiere GPU ≥20 GB VRAM).

Texto clínico de ejemplo:

    "The proband is a 4-year-old boy with renal failure and
    short stature. His parents are healthy and non-consanguineous."

Ejecución (requiere GPU):

    ~/.pga/envs/pga-core/bin/pga tool extract-hpo \
      --text "The proband is a 4-year-old boy with renal failure and short stature. His parents are healthy and non-consanguineous." \
      --out ./hpo \
      --backend local \
      --model mistral_24b

Salida esperada (con GPU):

    OK RDMA (local/mistral_24b): 3 términos HPO (de 5 entidades)
      entities_extracted: 5
      entities_verified: 3
      hpo_matched: 3
      hpo_unique: 3
      → ./hpo/hpo_terms.json

El fichero ./hpo/hpo_terms.json contiene:

    {
      "hpo_ids": ["HP:0000083", "HP:0004322", "HP:0012622"],
      "terms": [...]
    }

Sin GPU: usa los HPO de ejemplo manualmente. Define la variable
HPO una vez y úsala en los comandos siguientes:

    export HPO="HP:0000083,HP:0004322"

El resto del tutorial funciona igual.

## Paso 5. Rankear enfermedades candidatas

Con los términos HPO, LIRICAL prioriza enfermedades. Funciona sin GPU
(solo CPU) y tarda ~3 segundos.

    ~/.pga/envs/pga-core/bin/pga tool phenotype-ranking \
      --hpo "HP:0000083,HP:0004322" \
      --out ./rank \
      --top-n 10

Salida esperada:

    OK LIRICAL: 8525 enfermedades rankeadas, top-10 incluido en el report
      TSV: ./rank/lirical_<hash>.tsv
      Top 5:
          1  Alstrom syndrome  (OMIM:203800)  post=3,57%
          2  Galloway-Mowat syndrome 3  (OMIM:617729)  post=3,57%
          3  Cranioectodermal dysplasia  (OMIM:218330)  post=3,57%
          4  Muckle-Wells syndrome  (OMIM:191900)  post=3,57%
          5  Short-Rib thoracic dysplasia 10  (OMIM:615630)  post=3,57%

Las enfermedades relacionadas con NPHP1 (nefronoptisis) aparecen más
abajo en el ranking (~posición 57), porque los HPO de entrada son
genéricos. Ver docs/phenotype_layer.md para entender por qué.

El fichero ./rank/phenotype_ranking_report.json contiene:

    {
      "hpo_ids": ["HP:0000083", "HP:0004322"],
      "mode": "phenotype_only",
      "total_diseases_ranked": 8525,
      "top_n_included": 10
    }

## Paso 6. Cortar los BAMs al ROI de NPHP1

Para el análisis de variantes dirigido, cortamos solo el ROI de NPHP1.
El sistema no procesa el exoma completo.

    for s in c2 f1 m1; do
      ~/.pga/envs/pga-core/bin/pga tool fetch-roi \
        --bam cohort/$s.bam \
        --region chr2_roi:1-25000 \
        --out cohort/${s}_roi.bam
    done

Salida esperada:

    OK cohort/c2_roi.bam (NNN bytes, índice: c2_roi.bam.bai)
    OK cohort/f1_roi.bam (NNN bytes, índice: f1_roi.bam.bai)
    OK cohort/m1_roi.bam (NNN bytes, índice: m1_roi.bam.bai)

## Paso 7. QC de los sub-BAMs

Verificamos integridad, índice, orden y sample name.

    ~/.pga/envs/pga-core/bin/pga tool qc \
      --bam cohort/c2_roi.bam \
      --level counts \
      --expected-sample C2

Salida esperada:

    QC OK
    Header  SO=coordinate  refs=1  RG=1  SM=['C2']
    Counts  {"flagstat": {"mapped": 5000, ...}, "idxstats": {"chr2_roi": 5000}}

## Paso 8. Joint calling del trío

Llamamos variantes para cada individuo del trío contra GRCh38, y
fusionamos los VCFs con bcftools merge.

    for s in c2 f1 m1; do
      su=$(echo $s | tr '[:lower:]' '[:upper:]')
      ~/.pga/envs/pga-core/bin/pga tool call-variants \
        --bam cohort/${s}_roi.bam \
        --reference cohort/nphp1_ref.fa \
        --out cohort/${s}.vcf.gz \
        --region chr2_roi:1-25000 \
        --min-qual 0 --min-dp 1
    done

Salida esperada:

    OK cohort/c2.vcf.gz (26 variantes passing de 26, índice: c2.vcf.gz.tbi)
    OK cohort/f1.vcf.gz (29 variantes passing de 29, índice: f1.vcf.gz.tbi)
    OK cohort/m1.vcf.gz (23 variantes passing de 23, índice: m1.vcf.gz.tbi)

Renombramos el sample de cada VCF y fusionamos:

    for s in c2 f1 m1; do
      su=$(echo $s | tr '[:lower:]' '[:upper:]')
      ~/.pga/envs/pga-hts/bin/bcftools view cohort/${s}.vcf.gz | \
        sed "s/^\(#CHROM.*\)\tSAMPLE\$/\1\t$su/" | \
        ~/.pga/envs/pga-hts/bin/bcftools view -Oz -o cohort/${s}_named.vcf.gz
      ~/.pga/envs/pga-hts/bin/bcftools index -t cohort/${s}_named.vcf.gz
    done

    ~/.pga/envs/pga-hts/bin/bcftools merge -m all \
      cohort/c2_named.vcf.gz cohort/f1_named.vcf.gz cohort/m1_named.vcf.gz \
      -Oz -o cohort/trio.vcf.gz
    ~/.pga/envs/pga-hts/bin/bcftools index -t cohort/trio.vcf.gz

Verificación:

    ~/.pga/envs/pga-hts/bin/bcftools view -H cohort/trio.vcf.gz | wc -l

Salida esperada:

    31

El VCF del trío contiene 31 variantes (unión de las 26 + 29 + 23
individuales).

## Paso 9. Crear el .fam del trío

Para los filtros mendelianos necesitamos un .fam con la estructura
familiar. Formato: FID IID PID MID sex phenotype (1=sano, 2=afectado).

    cat > cohort/trio.fam << 'FAM'
    F1	F1	0	0	1	1
    M1	M1	0	0	2	1
    C2	C2	F1	M1	1	2
    FAM

Verificación:

    cat cohort/trio.fam

Salida esperada:

    F1	F1	0	0	1	1
    M1	M1	0	0	2	1
    C2	C2	F1	M1	1	2

F1 = padre varón sano, M1 = madre mujer sana, C2 = hijo varón afectado.

## Paso 10. Filtros mendelianos

Aplicamos los 4 modelos de segregación:

    ~/.pga/envs/pga-core/bin/pga tool mendelian-filter \
      --vcf cohort/trio.vcf.gz \
      --out cohort/mendel \
      --pedigree-fam cohort/trio.fam \
      --proband C2 \
      --no-filter-by-affected \
      --min-dp 5 --min-qual 0

Salida esperada:

    OK de_novo=0  auto_rec_hom=0  auto_dom=26  x_linked_rec=0
      de_novo: 0
      auto_rec_hom: 0
      auto_dom: 26
      x_linked_rec: 0

Interpretación clínica:

- **de_novo = 0**: ninguna variante presente solo en el hijo.
- **auto_rec_hom = 0**: ninguna variante homocigota recesiva.
- **auto_dom = 26**: 26 variantes donde el hijo y al menos un
  progenitor son portadores. Coherente con un patrón dominante
  o con penetrancia incompleta.
- **x_linked_rec = 0**: no hay cromosoma X en el fixture.

El reporte completo está en cohort/mendel/mendelian_report.json.

## Paso 11. (Opcional) Pangenómica local

El fixture contiene un cluster de SNVs MENA alrededor de la variante
de interés. Comparamos el calling lineal contra GRCh38 con el calling
contra un grafo local del ROI.

    ~/.pga/envs/pga-core/bin/pga tool build-local-graph \
      --reference cohort/nphp1_ref.fa \
      --cohort cohort/cohort.vcf.gz \
      --out cohort/roi.vg

    ~/.pga/envs/pga-hts/bin/samtools fastq -@ 2 cohort/c2_roi.bam > cohort/c2.fastq

    ~/.pga/envs/pga-core/bin/pga tool align-to-graph \
      --graph cohort/roi.vg \
      --reads cohort/c2.fastq \
      --out cohort/c2.gam \
      --pack cohort/c2.pack

    ~/.pga/envs/pga-core/bin/pga tool call-from-graph \
      --graph cohort/roi.vg \
      --xg cohort/roi.xg \
      --pack cohort/c2.pack \
      --out cohort/c2_graph.vcf.gz \
      --sample-name C2

Nota: vg call genera el VCF con sample "SAMPLE" por defecto. Usamos
--sample-name C2 para que el sample coincida con el del calling lineal
y compare-vcfs pueda emparejar las variantes.

Comparamos los dos callings:

    ~/.pga/envs/pga-core/bin/pga tool compare-vcfs \
      --baseline cohort/c2.vcf.gz \
      --candidate cohort/c2_graph.vcf.gz \
      --out cohort/delta.vcf.gz \
      --report cohort/delta.json \
      --sample C2

Salida esperada:

    OK /home/<user>/pga_quickstart/cohort/delta.vcf.gz
      consistent: 18
      discordant: 6
      recovered: 2
      lost: 2

Interpretación:

- **consistent (18)**: variantes con genotipo idéntico en ambos callings.
- **recovered (2)**: variantes que el calling contra el grafo detecta
  pero el calling lineal contra GRCh38 pierde. Candidatas a falsos
  negativos por reference bias.
- **lost (2)**: variantes detectadas por el lineal pero no por el grafo.
- **discordant (6)**: mismo locus, genotipo distinto entre los dos
  callings. Alelos con mapeo ambiguo a GRCh38.

El grafo recupera variantes que el calling lineal pierde, y absorbe
los alelos comunes MENA como parte de la referencia local. Ver
docs/demo_nphp1_mena.md para el análisis completo.

## Resumen

En 5 pasos principales has ejecutado el pipeline completo sobre el
fixture NPHP1:

1. Generado un trío sintético con frecuencias reales de gnomAD MID.
2. Extraído HPO del texto clínico (o usado HPO de ejemplo).
3. Rankeado enfermedades candidatas con LIRICAL.
4. Cortado los BAMs al ROI de NPHP1.
5. Aplicado filtros mendelianos al trío.
6. (Opcional) Comparado calling lineal vs grafo pangenómico local.

Los ficheros generados están en ~/pga_quickstart/cohort/ y
~/pga_quickstart/rank/.

## Estructura de los ficheros generados

    ~/pga_quickstart/
    ├── cohort/                    cohorte + trío
    │   ├── cohort.vcf.gz          cohort VCF (10 individuos)
    │   ├── pedigree.fam           estructura familiar
    │   ├── expected.json          variantes conocidas
    │   ├── nphp1_ref.fa           referencia del ROI
    │   ├── c2.bam, f1.bam, m1.bam BAMs del trío
    │   ├── c2_roi.bam, ...        sub-BAMs del ROI
    │   ├── c2.vcf.gz, f1.vcf.gz, m1.vcf.gz   VCFs individuales
    │   ├── trio.vcf.gz            VCF conjunto del trío
    │   ├── trio.fam               .fam del trío
    │   ├── mendel/                resultados de mendelian_filter
    │   │   ├── de_novo.vcf.gz
    │   │   ├── auto_rec_hom.vcf.gz
    │   │   ├── auto_dom.vcf.gz
    │   │   ├── x_linked_rec.vcf.gz
    │   │   └── mendelian_report.json
    │   ├── roi.vg, roi.xg         grafo pangenómico local
    │   ├── c2.gam, c2.pack        alineamiento contra el grafo
    │   ├── c2_graph.vcf.gz        calling contra el grafo
    │   └── delta.vcf.gz           delta lineal vs grafo
    ├── hpo/                       extracción HPO (si hay GPU)
    │   └── hpo_terms.json
    └── rank/                      ranking de enfermedades
        ├── lirical_<hash>.tsv
        └── phenotype_ranking_report.json

## Qué hacer a continuación

- **Analizar tus propios BAMs**: sustituye cohort/*.bam por tus
  ficheros, y `cohort/cohort.vcf.gz` por tu cohort de referencia
  (gnomAD MID, Al Mena, KFSHRC, o cualquier VCF anotado).
- **Explorar la documentación**: docs/ tiene el diseño conceptual
  completo, desde la virtualización hasta la capa agéntica.
- **Ejecutar los tests**: `pytest tests/ -v` en el entorno pga-core
  valida los 26 tests de integración.
- **Leer las demos reproducibles**: docs/demo_nphp1*.md contiene
  los 4 experimentos científicos que validan el sistema contra
  escenarios sintéticos y cohortes MENA.

## Recursos

- **README** — instalación y uso rápido.
- **docs/installation.md** — instalación detallada y troubleshooting.
- **docs/use_cases.md** — cuatro escenarios de uso (clínico,
  fenotipo, poblacional, GWAS).
- **docs/phenotype_layer.md** — diseño de la capa de fenotipo
  (RDMA, LIRICAL, requisitos de hardware).
- **docs/local_pangenome.md** — diseño de la pangenómica local.
- **docs/mendelian_segregation.md** — filtros mendelianos.
- **docs/decision_not_pipeline.md** — por qué este sistema no es
  un pipeline más.

## Limitaciones conocidas

- **La extracción HPO (RDMA) requiere GPU ≥20 GB VRAM.** Sin GPU,
  usa HPO manuales. Ver docs/phenotype_layer.md sección 14.
- **LIRICAL requiere Java 17+.** Ver docs/lirical_setup.md.
- **El fixture es sintético.** Los BAMs y los genotipos del trío
  son simulados a partir de frecuencias reales. Para uso clínico
  real, sustituye por BAMs y cohorts reales.
- **No alineamos.** PublicGenomicAgent empieza con BAMs ya
  alineados a hg38.
