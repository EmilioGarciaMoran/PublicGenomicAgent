# PublicGenomicAgent

Agente genómico **clínico** conversacional, quirúrgico y virtualizado.

Dado un fenotipo, una familia y BAMs ya alineados, analiza solo las
regiones relevantes (ROIs) y devuelve un VCF delta con variantes
candidatas priorizadas. Optimizado para hardware modesto: portátiles,
discos externos, sin clúster.

## Qué es esto

PublicGenomicAgent no es un pipeline. Un pipeline produce ficheros;
este sistema produce decisiones.

Entrada: un BAM ya alineado + una consulta clínica (fenotipo,
familia, sospecha diagnóstica) + opcionalmente un texto clínico.

Salida: una decisión cuantificada sobre variantes candidatas,
con estadística mendeliana, pangenomización local, y trazabilidad
completa del proceso.

El VCF es el lenguaje intermedio. La decisión es el producto.

## Enfoque: clínico, no genómico

| Aspecto | Pipeline genómico | PublicGenomicAgent |
|---|---|---|
| Punto de partida | FASTQ / exoma completo | BAM + fenotipo |
| Scope | Todo el genoma | ROI dirigido por clínica |
| Tiempo por caso | Horas | Segundos por ROI |
| Salida | VCF exhaustivo | VCF delta priorizado |
| Usuario | Bioinformático | Genetista clínico |
| Validación | Sensibilidad global | Utilidad diagnóstica |

## No alineamos

El sistema empieza con **BAMs ya alineados** a hg38, producidos por
vuestro pipeline habitual (DRAGEN, BWA+GATK, minimap2, etc.).
Asumimos que el trabajo pesado de alineamiento ya está hecho, en
servidor o en cloud.

El valor del sistema está en la capa clínica: interpretación dirigida
por fenotipo sobre BAMs existentes, con pangenomización local para
corregir el reference bias poblacional.

## Pipeline

Texto clínico
    ↓
extract_hpo (RDMA, GPU)  →  HPO terms
    ↓
phenotype_ranking (LIRICAL, CPU)  →  enfermedades rankeadas
    ↓
fetch_roi → qc → call_variants → local_pangenome → build_local_graph
    ↓
mendelian_filter → plink_validate
    ↓
compare_vcfs → delta

Cada paso es una tool del CLI. Todas son deterministas, auditables
y se pueden encadenar.

## Instalación

Requisitos:
- Linux o macOS
- micromamba en PATH
- Java 17+ (para LIRICAL)
- Opcional: GPU NVIDIA ≥20 GB (para extracción HPO con RDMA)

### Bootstrap desde cero

    git clone https://github.com/EmilioGarciaMoran/PublicGenomicAgent
    cd PublicGenomicAgent

    # 1. Crear el entorno base (una sola vez)
    micromamba create -y -p ~/.pga/envs/pga-core -f envs/pga-core.yml
    micromamba run -p ~/.pga/envs/pga-core pip install -e .

    # 2. Crear los entornos de herramientas
    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-hts
    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-pangenome
    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-mendelian
    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-phenotype

    # 3. Bootstrap de la capa de fenotipo (datos de RDMA + LIRICAL)
    ~/.pga/envs/pga-core/bin/pga phenotype bootstrap

    # 4. Verificar
    ~/.pga/envs/pga-core/bin/pga env verify

Con eso el sistema está operativo. Ver docs/installation.md para
detalles y solución de problemas.

## Uso rápido

### Capa genómica (sin GPU)

    pga tool fetch-roi --bam patient.bam --region chr2:110000000-110025000 --out roi.bam
    pga tool qc --bam roi.bam --level counts --expected-sample PAT001
    pga tool call-variants --bam roi.bam --reference hg38_roi.fa --out patient.vcf.gz
    pga tool mendelian-filter --vcf trio.vcf.gz --pedigree-fam trio.fam --proband C2

### Capa de fenotipo

    # Extracción de HPO desde texto (requiere GPU ≥20 GB)
    pga tool extract-hpo --text "patient with renal failure" --out ./hpo

    # Priorización de enfermedades desde HPO (solo CPU)
    pga tool phenotype-ranking --hpo "HP:0000083,HP:0004322" --out ./rank

### Pangenomización local

    pga tool local-pangenome --cohort cohort.vcf.gz --reference hg38_roi.fa \
      --region chr2_roi:1-25000 --out ./enriched --min-af 0.05

    pga tool build-local-graph --reference hg38_roi.fa \
      --cohort cohort.vcf.gz --out roi.vg

    pga tool compare-vcfs --baseline patient_grch38.vcf.gz \
      --candidate patient_pga.vcf.gz --out delta.vcf.gz --report delta.json

Ver docs/quickstart.md para un ejemplo end-to-end con el fixture
NPHP1 (trío con fallo renal y enanismo).

## Estructura

    envs/         manifiestos y registry de entornos
    src/          código del paquete publicgenomicagent
    docs/         documentación de diseño
    tests/        tests unitarios e integración
    scripts/      scripts de bootstrap y utilidades

## Documentación

Diseño conceptual:

- docs/virtualization.md — capa de entornos aislados.
- docs/local_pangenome.md — núcleo: pangenomización local.
- docs/decision_not_pipeline.md — por qué es distinto.
- docs/mendelian_segregation.md — la genética clínica como lógica.
- docs/phenotype_layer.md — capa de fenotipo.
- docs/audit_and_training.md — sesiones registradas y formación.

Guías de uso:

- docs/installation.md — instalación detallada y troubleshooting.
- docs/quickstart.md — ejemplo end-to-end.
- docs/use_cases.md — cuatro escenarios de uso.

Demos reproducibles:

- docs/demo_nphp1.md — experimento lineal.
- docs/demo_nphp1_graph.md — comparación lineal vs grafo.
- docs/demo_nphp1_mena.md — grafo recupera SNVs MENA.
- docs/demo_nphp1_trio.md — joint calling del trío.

Setup de herramientas externas:

- docs/lirical_setup.md — instalación y bug conocido.

## Estado

En desarrollo activo. 12 tools funcionales, 5 entornos aislados,
26 tests verdes.

Capa genómica (sin GPU): completa.

Capa de fenotipo (extracción HPO): requiere GPU ≥20 GB o servidor
dedicado. Ver docs/phenotype_layer.md sección 14.

## Licencia

MIT. Ver LICENSE.
