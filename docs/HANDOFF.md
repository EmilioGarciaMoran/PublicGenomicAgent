# Handoff — PublicGenomicAgent

**Última actualización**: 2026-10-03
**Estado**: Fase 2 al 100% (salvo item 2.2 bloqueado por descargas).

## Estado del repo

- **Último commit**: b6b3aa6
- **Tests**: 119 passed (85 unit, 34 integration)
- **CI**: verde (GitHub Actions)
- **Branch**: main, sync con origin/main

## Qué se ha hecho

### Sesión 2026-10-03 (última)
- Añadido `prioritize_variants` en el pipeline (paso 13 con --clinvar).
- Test end-to-end del pangenoma (`test_pangenome_pipeline.py`).
- Actualizados `docs/quickstart.md`, `docs/backlog.md`, README.

### Sesión 2026-10-02
- ... (resumen)

## Qué falta

### Fase 2 (pendiente)
- **2.2** Cohortes locales (gnomAD MENA, KSA) — bloqueado por descargas.

### Fase 3 (próxima)
- **3.1** Trío real (1000G, KSA001).
- **3.2** Benchmark lineal vs pangenoma (F1-score).
- **3.3** Pitch para KFSH&RC.
- **3.4** Quickstart con caso real.

### Fase 4 (futuro)
- ...

## Cómo retomar

```bash
cd ~/PublicGenomicAgent
git pull
/home/egarmo/bin/micromamba run -p ~/.pga/envs/pga-core pytest tests/ -q

# Pipeline completo (13 pasos)
/home/egarmo/bin/micromamba run -p ~/.pga/envs/pga-core \
  pga run-trio ... --clinvar ~/clinvar.vcf.gz
