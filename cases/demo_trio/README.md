# Caso demo: trío sintético NPHP1

Trío sintético (father / mother / proband) con una región de interés
en `chr2` que engloba el locus NPHP1. Los BAMs fueron generados con
`BAMMaker` (proyecto local separado) usando SPDI sobre un cromosoma
sintético de 200 kb.

## Contenido

- `manifest.json` — `CaseManifest` que consume `pga chat --case`.
- `bams/` — tres BAMs con sus índices `.bai`.
- `ref/` — `test_ref.fa` (200 kb, chr2) con su `.fai`.

## Uso

```bash
# Verificar que está todo
python scripts/make_demo_case.py

# Arrancar el chat sobre el caso
pga chat --case cases/demo_trio/manifest.json

# Auditar la sesión resultante
pga audit sessions/$(ls -t sessions/ | head -n 1)/session.json --verify
