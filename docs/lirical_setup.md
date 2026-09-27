# Setup de LIRICAL

Guía de instalación y configuración de LIRICAL en PublicGenomicAgent.

## Requisitos

- Java 17+ (probado con Java 21)
- ~400 MB de disco para el JAR + datos

## Instalación

    mkdir -p ~/.pga/cache/lirical
    cd ~/.pga/cache/lirical
    curl -L -O https://github.com/TheJacksonLaboratory/LIRICAL/releases/latest/download/lirical-cli-2.4.1-distribution.zip
    unzip -q lirical-cli-2.4.1-distribution.zip
    cd lirical-cli-2.4.1
    java -jar lirical-cli-2.4.1.jar download -d data

## Bug conocido: duplicados en datos HGNC (2026)

LIRICAL v2.4.1 falla al arrancar con los datos actuales de HGNC y NCBI
debido a dos duplicados:

    1. mim2gene_medgen: mismo NCBIGene ID para dos genes distintos
    2. hgnc_complete_set.txt: mismo Ensembl ID para dos genes distintos
       (ejemplo: ENSG00000230417 asignado a LINC00595 y LINC00856)

Síntoma:

    java.lang.IllegalStateException: Duplicate key NCBIGene:ENSG...

Causa: la librería phenol que usa LIRICAL construye mapas con
Collectors.toMap, que falla si hay claves duplicadas.

Fix: filtrar duplicados de ambos ficheros antes de ejecutar LIRICAL.

    cd ~/.pga/cache/lirical/lirical-cli-2.4.1

    # Filtrar mim2gene_medgen por NCBIGene ID (columna 2)
    awk -F'\t' '
      /^#/ { print; next }
      $2 == "-" { print; next }
      !seen[$2]++ { print }
    ' data/mim2gene_medgen > data/mim2gene_medgen.fixed
    mv data/mim2gene_medgen.fixed data/mim2gene_medgen

    # Filtrar hgnc_complete_set.txt por Ensembl ID (columna 19)
    awk -F'\t' '
      NR==1 { print; next }
      $19 == "" { print; next }
      !seen[$19]++ { print }
    ' data/hgnc_complete_set.txt > data/hgnc_complete_set.fixed
    mv data/hgnc_complete_set.fixed data/hgnc_complete_set.txt

Este bug debería reportarse en el repositorio de LIRICAL.

## Verificación

    cd ~/.pga/cache/lirical/lirical-cli-2.4.1

    java -jar lirical-cli-2.4.1.jar prioritize \
      -p "HP:0000083,HP:0004322" \
      -d ./data \
      -x test \
      -o /tmp/lirical_test \
      -f tsv

Esperado: análisis completado en <5 segundos, con un TSV con
~8.500 enfermedades rankeadas.

## Uso desde PublicGenomicAgent

    ~/.pga/envs/pga-core/bin/pga tool phenotype-ranking \
      --hpo "HP:0000083,HP:0004322" \
      --out /tmp/phenotype_ranking_test \
      --top-n 20

## Modo genotipo-aware

    ~/.pga/envs/pga-core/bin/pga tool phenotype-ranking \
      --hpo "HP:0000083,HP:0004322" \
      --vcf patient.vcf.gz \
      --assembly hg38 \
      --out /tmp/phenotype_ranking_ga

Requiere haber configurado las bases de datos de Exomiser (hg19/hg38).
