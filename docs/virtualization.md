# Virtualización de PublicGenomicAgent

Cómo PublicGenomicAgent aísla sus herramientas bioinformáticas para evitar
el caos de versiones, las instalaciones zombies y las dependencias globales.

## Principios

1. Nada global. Ni pip global, ni conda base, ni sudo apt install.
   Todo binario vive en ~/.pga/envs/.
2. Un entorno por familia de herramientas. samtools no comparte prefix
   con vep ni con mafft.
3. Manifiestos versionables, entornos desechables. La fuente de verdad
   es envs/*.yml.
4. Hash de manifiesto como identificador de estado. Un entorno está
   actualizado solo si su lock coincide con el hash del manifiesto.
5. Ejecución por ruta absoluta. Los binarios se invocan directamente.
6. El paquete Python propio vive en un solo entorno. publicgenomicagent
   se instala en pga-core vía installs_self: true.

## Familias de entornos

| Entorno     | Contenido                                  | Estado      |
|-------------|--------------------------------------------|-------------|
| pga-core    | Python, CLI, agente, estado, knowledge     | activo      |
| pga-hts     | htslib (samtools, bcftools, tabix, bgzip)  | activo      |
| pga-annot   | VEP, snpEff, pysam, cyvcf2                 | planificado |
| pga-phylo   | MAFFT, FastTree, IQ-TREE, biopython        | planificado |
| pga-report  | igv-reports, jinja2                        | planificado |

## Layout en disco

Estado de runtime bajo ~/.pga/ (o $PGA_ROOT si está definido):

    ~/.pga/
      envs/     prefijos micromamba
      cache/    referencias públicas (hg38, HPO, ClinVar)
      locks/    ficheros <env>.lock.json
      logs/     logs de bootstrap

En el repositorio solo viven manifiestos, registry y el código de
src/publicgenomicagent/env/.

## Manifiestos

Cada manifiesto declara dependencias exactas de un entorno:

    name: pga-hts
    channels:
      - conda-forge
      - bioconda
    dependencies:
      - htslib=1.19.1
      - samtools=1.19.2
      - bcftools=1.19

Reglas:

- Versiones pinneadas con =X.Y.Z siempre que sea posible.
- Canales en orden: conda-forge primero, bioconda después.
- Nada de -e . ni instalaciones del paquete propio en el manifiesto.
- El nombre debe coincidir con la clave en registry.yaml.

## Registry

envs/registry.yaml es el índice único del sistema. Declara entornos y
herramientas:

    envs:
      pga-core:
        manifest: envs/pga-core.yml
        description: "Agente, CLI, orquestación, estado"
        installs_self: true
      pga-hts:
        manifest: envs/pga-hts.yml
        description: "htslib (samtools, bcftools, tabix, bgzip)"

    tools:
      samtools:
        env: pga-hts
        binary: samtools
        version_flag: "--version"
        min_version: "1.19"

Campos:

- installs_self: si true, tras crear el entorno se instala
  publicgenomicagent en editable. Solo debe estar activo en pga-core.
- tools.<name>.env: entorno donde vive el binario.
- tools.<name>.version_flag: cómo pedir la versión.
- tools.<name>.min_version: umbral mínimo (documental por ahora).

## Ciclo de vida del bootstrap

pga env bootstrap <nombre> decide entre cuatro caminos:

| Estado                                | Acción                                  |
|---------------------------------------|-----------------------------------------|
| No existe el prefix                   | Crear entorno, instalar pkg si aplica   |
| Existe prefix, no hay lock            | Registrar lock, instalar pkg si aplica  |
| Existe lock, hash distinto            | Recrear entorno, instalar pkg si aplica |
| Existe lock, hash igual               | No hacer nada (idempotente)             |

## Resolución de herramientas

Cuando el agente pide samtools, ToolRuntime:

1. Consulta registry.yaml y ve que samtools vive en pga-hts.
2. Verifica que pga-hts está creado.
3. Construye la ruta ~/.pga/envs/pga-hts/bin/samtools.
4. Ejecuta por ruta absoluta con subprocess.run.

## Verificación

pga env verify recorre entornos y herramientas y reporta:

- ok — el binario responde y su versión es coherente.
- stale_manifest — el manifiesto cambió.
- missing_env — el entorno declarado no existe.
- missing_binary — el entorno existe pero el binario no aparece.

## Quirks conocidos

### micromamba run intercepta --version

En micromamba 2.x, micromamba run -p <prefix> -- python --version devuelve
la versión de micromamba en lugar de la del binario. Por eso el runtime
ejecuta binarios por ruta absoluta, sin pasar por micromamba run.

### Algunos binarios imprimen --version en stderr

samtools y bcftools escriben el banner de versión en stderr. Por eso
ToolRuntime.capture() unifica stdout + stderr antes de parsear.

### tabix y bgzip no son paquetes sueltos

No existen como paquetes independientes en bioconda con versión propia.
Se instalan como parte de htslib. No añadirlos al manifiesto como
dependencias separadas.

## Reproducibilidad

Cada sesión de análisis emitirá un session.json con:

- Hash de cada manifiesto usado.
- Versión de cada binario invocado.
- Hash del genoma de referencia.
- Hash de la base de conocimiento (DuckDB).
- Hash de los ficheros de entrada.
- Parámetros del agente (prompt, temperatura, modelo LLM).

Con eso, cualquier sesión puede reconstruirse bit a bit sobre un clon
limpio del repositorio y los mismos datos de entrada.
