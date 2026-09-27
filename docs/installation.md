# Instalación detallada

Guía completa de instalación de PublicGenomicAgent, con explicación de
qué hace cada paso y cómo resolver los problemas más comunes.

Para el bootstrap rápido, ver el README.

## Requisitos

Obligatorios:

- Linux o macOS
- micromamba en PATH (o MICROMAMBA_BIN definido)
- Java 17+ (LIRICAL lo requiere)

Opcionales pero recomendados:

- GPU NVIDIA con ≥20 GB VRAM (para extracción HPO con RDMA)
- 10 GB de disco libre en ~/.pga/ (entornos + datos + cache)

## Qué hace la instalación

PublicGenomicAgent se instala en tres capas independientes:

1. **Entornos aislados** (micromamba): 5 entornos con herramientas
   bioinformáticas. No tocan el Python del sistema.
2. **Paquete Python propio** (publicgenomicagent): instalado en el
   entorno pga-core en modo editable.
3. **Datos de referencia** (cache): vector stores, LIRICAL, HPO,
   cohortes. Se descargan a ~/.pga/cache/.

Cada capa se puede actualizar o borrar sin afectar a las demás.

## Instalación paso a paso

### Paso 1. Clonar el repositorio

    git clone https://github.com/EmilioGarciaMoran/PublicGenomicAgent
    cd PublicGenomicAgent

### Paso 2. Crear el entorno base (pga-core)

Este es el único paso manual. A partir de aquí, todo se hace con
el CLI pga.

    micromamba create -y -p ~/.pga/envs/pga-core -f envs/pga-core.yml
    micromamba run -p ~/.pga/envs/pga-core pip install -e .

El entorno pga-core contiene Python 3.11, typer, rich, pydantic,
pyyaml y el propio paquete publicgenomicagent en modo editable.

### Paso 3. Crear los entornos de herramientas

    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-hts
    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-pangenome
    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-mendelian
    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-phenotype

Cada bootstrap descarga e instala las dependencias del entorno
correspondiente en ~/.pga/envs/<nombre>/. Son ~500 MB por entorno.

### Paso 4. Bootstrap de la capa de fenotipo

    ~/.pga/envs/pga-core/bin/pga phenotype bootstrap

Descarga:

- Vector stores de RDMA (~300 MB) desde GitHub Releases
- LIRICAL 2.4.1 + datos (~450 MB) con el fix de duplicados aplicado

### Paso 5. Verificar

    ~/.pga/envs/pga-core/bin/pga env list
    ~/.pga/envs/pga-core/bin/pga env verify

Esperado: 5 entornos declarados, 5 herramientas verificadas.

## Estructura en disco

Todo lo que instala PublicGenomicAgent vive bajo ~/.pga/:

    ~/.pga/
    ├── envs/              entornos micromamba (5, ~500 MB cada uno)
    │   ├── pga-core/
    │   ├── pga-hts/
    │   ├── pga-pangenome/
    │   ├── pga-mendelian/
    │   └── pga-phenotype/
    ├── cache/             datos descargables
    │   ├── reference/     hg38 y secuencias UCSC
    │   ├── gnomad/        cohortes gnomAD cacheadas
    │   ├── rdma/          vector stores de RDMA (~300 MB)
    │   ├── lirical/       LIRICAL + datos (~450 MB)
    │   └── models/        modelos LLM (Mistral 24B, opcional)
    ├── locks/             hashes de manifiestos por entorno
    └── logs/              trazas de bootstrap

Nada de esto entra en el repositorio. Cada máquina que clona el repo
construye su propio ~/.pga/ desde cero.

## Verificación

### Ver los entornos declarados

    ~/.pga/envs/pga-core/bin/pga env list

Salida esperada:

    pga-core       creado
    pga-hts        creado
    pga-pangenome  creado
    pga-mendelian  creado
    pga-phenotype  creado

### Verificar herramientas

    ~/.pga/envs/pga-core/bin/pga env verify

Salida esperada:

    pga-core       python     ok    3.11.16
    pga-hts        samtools   ok    1.19.2
    pga-hts        bcftools   ok    1.19
    pga-pangenome  vg         ok    1.76.1
    pga-mendelian  plink      ok    1.90

### Verificar la capa de fenotipo

    ~/.pga/envs/pga-phenotype/bin/python -c "from rdma.hpo.extractor import PhenotypeExtractor; print('RDMA OK')"

    ls ~/.pga/cache/rdma/vector_stores/

Salida esperada:

    RDMA OK
    G2GHPO_metadata_medembed.npy
    rd_orpha_medembed.npy

## Solución de problemas

### micromamba no encontrado

Síntoma:

    RuntimeError: micromamba no encontrado

Solución: instala micromamba o define la variable de entorno:

    export MICROMAMBA_BIN=/ruta/a/micromamba

o edita ~/.bashrc para que persista.

### Java no instalado

Síntoma:

    RuntimeError: Java no encontrado. LIRICAL requiere Java 17+.

Solución:

    sudo apt install openjdk-17-jre        # Ubuntu / Debian
    brew install openjdk@17                 # macOS

Verifica con:

    java -version

### GPU no disponible para RDMA

Síntoma:

    RuntimeError: Backend 'local' requiere GPU NVIDIA con CUDA.
    No se detectó CUDA disponible.

Esto es esperado si tu máquina no tiene GPU NVIDIA. La capa genómica
funciona sin GPU. La extracción HPO (extract_hpo) requiere GPU ≥20 GB.

Alternativas:

- Ejecutar extract_hpo en un servidor con GPU.
- Usar el backend openrouter o api (solo para desarrollo, no para
  datos clínicos reales).

Ver docs/phenotype_layer.md sección 14 para detalles.

### Error de duplicados en LIRICAL

Síntoma:

    java.lang.IllegalStateException: Duplicate key NCBIGene:ENSG00000230417

Causa: LIRICAL 2.4.1 con datos HGNC/NCBI de 2026 tiene duplicados en
mim2gene_medgen y hgnc_complete_set.txt.

Solución: el bootstrap lo arregla automáticamente. Si estás usando
LIRICAL manualmente, ver docs/lirical_setup.md para el fix.

### Espacio en disco insuficiente

Síntoma:

    OSError: [Errno 28] No space left on device

Causa: los 5 entornos + cache ocupan ~10 GB. Verifica con:

    du -sh ~/.pga/

Solución: libera espacio o mueve ~/.pga/ a otro disco:

    export PGA_ROOT=/ruta/a/otro/disco/.pga

### Entorno corrupto tras interrupción

Si un bootstrap se interrumpe a mitad (Ctrl+C, caída de red), el
entorno puede quedar en estado inconsistente.

Solución: forzar reconstrucción:

    ~/.pga/envs/pga-core/bin/pga env bootstrap <entorno> --force

o borrar el entorno y relanzar:

    rm -rf ~/.pga/envs/<entorno> ~/.pga/locks/<entorno>.lock.json
    ~/.pga/envs/pga-core/bin/pga env bootstrap <entorno>

### Los tests fallan tras actualizar el repo

Si has hecho git pull y los tests fallan, probablemente algún
manifiesto ha cambiado.

Solución: reconstruir los entornos afectados.

    ~/.pga/envs/pga-core/bin/pga env verify

Verás qué entornos tienen stale_manifest. Reconstruye cada uno:

    ~/.pga/envs/pga-core/bin/pga env bootstrap <entorno>

## Desinstalación

Toda la instalación vive bajo ~/.pga/. Para borrarla:

    rm -rf ~/.pga/

Eso elimina los entornos, la cache, los locks y los logs. No toca
el repositorio ni el Python del sistema.

Para borrar también el repo:

    rm -rf ~/PublicGenomicAgent
