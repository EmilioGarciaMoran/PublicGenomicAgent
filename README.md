# PublicGenomicAgent

Agente genómico **clínico** conversacional, quirúrgico y virtualizado.

Dado un fenotipo, una familia y BAMs ya alineados, analiza solo las
regiones relevantes (ROIs) y devuelve un VCF delta con variantes
candidatas priorizadas. Optimizado para hardware modesto: portátiles,
discos externos, sin clúster.

## Enfoque: clínico, no genómico

PublicGenomicAgent no es un pipeline genómico más. Es una herramienta
clínica que usa genómica como sustrato.

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

## Núcleo: pangenomización local

En poblaciones infrarrepresentadas (Middle East, minorías), GRCh38
produce falsos negativos por reference bias. El sistema construye
un **consenso enriquecido del ROI** que incorpora alelos comunes de
una cohorte poblacional, y llama variantes contra ese consenso en
lugar de contra GRCh38.

La métrica de éxito no es "el consenso es bonito". Es el **VCF
diferencial**: `patient_grch38.vcf` vs `patient_enriched.vcf`. Las
variantes que aparecen solo en el segundo son candidatas a falsos
negativos recuperados.

Ver `docs/local_pangenome.md` para el diseño completo.

## Estado

En desarrollo activo. Implementado:

- Capa de virtualización con micromamba (`pga env`).
- `fetch_roi`: corta BAMs por rango genómico.
- `qc_bam`: control de calidad estructural y de conteos.
- `call_variants`: variant calling quirúrgico (single-sample).
- `local_pangenome`: consenso enriquecido del ROI.
- `compare_vcfs`: benchmarking VCF vs VCF (delta categorizado).

Pendiente:

- `mendelian`: filtros según pedigree.
- Phasing por trío.
- Simulador de BAM del paciente.
- Capa agéntica (LLM) y UI conversacional.
- Cohortes locales reales (KFSHRC, Al Mena).

## Requisitos

- Linux o macOS
- [micromamba](https://mamba.readthedocs.io/) en `PATH`, o
  `MICROMAMBA_BIN=/ruta/a/micromamba`
- Espacio en disco para entornos (~500 MB por familia)
- BAMs ya alineados a hg38 (no incluidos)

## Arranque desde cero

Crear el entorno base (único paso manual):

    micromamba create -y -p ~/.pga/envs/pga-core -f envs/pga-core.yml
    micromamba run -p ~/.pga/envs/pga-core pip install -e .

A partir de ahí, todo vía CLI:

    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-hts
    ~/.pga/envs/pga-core/bin/pga env list
    ~/.pga/envs/pga-core/bin/pga env verify

## Uso típico

    # Cortar el BAM del paciente al ROI clínico
    pga tool fetch-roi --bam patient.bam \
      --region chr2:110000000-110025000 --out roi.bam

    # Control de calidad
    pga tool qc --bam roi.bam --level counts --expected-sample PAT001

    # Construir la referencia enriquecida del ROI
    pga tool local-pangenome --cohort gnomad_mid.vcf.gz \
      --reference hg38_roi.fa --region chr2_roi:1-25000 \
      --out enriched/ --min-af 0.05 --af-field AF_MID

    # Llamar variantes contra GRCh38 (baseline) y contra el enriquecido
    pga tool call-variants --bam roi.bam --reference hg38_roi.fa \
      --out patient_grch38.vcf.gz
    pga tool call-variants --bam roi.bam --reference enriched/enriched_roi.fa \
      --out patient_enriched.vcf.gz

    # Comparar los dos callings
    pga tool compare-vcfs --baseline patient_grch38.vcf.gz \
      --candidate patient_enriched.vcf.gz \
      --out delta.vcf.gz --report delta.json

## Estructura

    envs/         manifiestos y registry de entornos
    src/          código del paquete publicgenomicagent
    docs/         documentación de diseño
    tests/        tests unitarios e integración

## Documentación

- `docs/virtualization.md` — capa de entornos aislados.
- `docs/local_pangenome.md` — núcleo: pangenomización local.
- `docs/agentic_layer.md` — restricción de diseño para el agente LLM.

## Cita

Si usas este software en investigación, ver `CITATION.cff` (pendiente).
