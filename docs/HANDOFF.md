# Handoff — PublicGenomicAgent

**Última actualización**: 2026-10-05
**Estado**: Fase 2 al 100% (salvo item 2.2 bloqueado por descargas).

## Estado del repo

- **Último commit**: 0c31431 (benchmark: pangenome wins on synthetic SD)
- **Tests**: 120 passed (85 unit, 35 integration)
- **CI**: verde (GitHub Actions)
- **Branch**: main, sync con origin/main

## Qué se ha hecho

### Sesión 2026-10-05 (última)
- **Benchmark en duplicación segmentaria sintética** (chr1, tandem dup de 100 bp):
  * Lineal (default): 0 variantes detectadas.
  * Lineal (sin filtro MAPQ): 0 variantes detectadas.
  * Pangenoma (vg): **1 variante, GT=1/1, DP=49, GQ=132**.
  * `docs/benchmark_duplications.md`.
- **Limpieza de disco**: de 18 GB a 30 GB libres (borrados BAMs de Osteopetrosis,
  artefactos intermedios del benchmark, `pga-phenotype`).
- **Test de integración del benchmark** lineal vs pangenoma.

### Sesión 2026-10-04
- Fix del pipeline con ClinVar: `_rule_annotate` corre antes de `_rule_mendelian`.
- `_rule_mendelian` ahora prefiere `state.vcfs['annotated']`.
- `_rule_prioritize` ahora prefiere `auto_rec_hom.vcf.gz`.
- Primer run real: Osteopetrosis trio (chr8, hg19) en 19.6s con `--clinvar`.
- `docs/demo_osteo_real.md` documenta el hito.
- Primer benchmark lineal vs pangenoma (102 vs 91 variantes, 76 compartidas).
- `docs/benchmark_linear_vs_pangenome.md`.

### Sesión 2026-10-03
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
