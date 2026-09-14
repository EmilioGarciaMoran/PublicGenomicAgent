# PublicGenomicAgent

Agente genómico conversacional, quirúrgico y virtualizado.

Analiza datos genómicos clínicos dirigido por fenotipo, familia o cohorte,
extrayendo solo las regiones de interés (ROI) de los BAMs. Optimizado para
hardware modesto: portátiles, discos externos, sin clúster.

## Estado

En desarrollo. La capa de virtualización (`pga env ...`) es funcional.

## Arranque desde cero

Crear el entorno base (único paso manual):

    micromamba create -y -p ~/.pga/envs/pga-core -f envs/pga-core.yml
    micromamba run -p ~/.pga/envs/pga-core pip install -e .

A partir de ahí, todo vía CLI:

    ~/.pga/envs/pga-core/bin/pga env bootstrap pga-hts
    ~/.pga/envs/pga-core/bin/pga env list
    ~/.pga/envs/pga-core/bin/pga env verify

## Documentación

- docs/virtualization.md
