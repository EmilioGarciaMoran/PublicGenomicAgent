# Backlog

Work in progress and pending items, organised by priority and
impact. This file is the entry point for anyone joining the
project (including future-me).

Legend: ✅ done · ⚠️ partial · ❌ not started · 🚧 in progress

---

## Current state (summary)

| Layer | Status |
|---|---|
| Virtualization (micromamba, 5 envs) | ✅ |
| Tools (13 registered, Pydantic I/O) | ✅ |
| Agent state, loop, planners | ✅ |
| CLI (`env`, `tool`, `llm`, `plan`, `run-trio`) | ✅ |
| Pipeline end-to-end on synthetic trio (11 steps) | ✅ |
| Pipeline with --clinvar (12 steps) | ✅ |
| HTML report (with variant table) | ✅ |
| IGV session + control URLs | ✅ |
| Tests (77 unit + 30 integration) | ✅ |
| Privacy policy (`docs/privacy.md`) | ✅ |
| README (English, professional) | ✅ |
| License (Apache 2.0) | ✅ |
| GitHub public repo, sync'd | ✅ |
| Annotation layer (VEP, ClinVar, gnomAD) | ❌ |
| NL → CaseManifest resolution | ❌ |
| HPO extraction without GPU | ❌ |
| Clinical report (HTML/PDF) | ❌ |
| IGV integration | ❌ |
| Pangenome end-to-end test | ❌ |
| CI (GitHub Actions) | ❌ |
| `setup.sh` unified installer | ❌ |
| Benchmark vs linear / commercial | ❌ |

---

## Phase 1 — Consolidate (this week)

Close the operational gaps that make the project
self-contained and installable by a third party.

| # | Item | Effort | Impact | Status |
|---|---|---|---|---|
| 1.1 | `docs/backlog.md` (this file) | 15 min | Medium | 🚧 |
| 1.2 | `setup.sh` unified installer | 30 min | High | ❌ |
| 1.3 | GitHub Actions CI (unit + integration) | 30 min | Medium | ❌ |
| 1.4 | `docs/quickstart.md` update | 30 min | Medium | ❌ |
| 1.5 | `docs/known_issues.md` confirm BAQ fix | 15 min | Low | ❌ |

---

## Phase 2 — Clinical value (2-3 weeks)

The layer that turns the project from a technical
framework into a tool a geneticist would actually use.

| # | Item | Effort | Impact | Status |
|---|---|---|---|---|
| 2.1 | `annotate_variants` tool (ClinVar) | 3-4h | High | ✅ |
| 2.1b | `--clinvar` auto-annotation in `run-trio` | 30min | High | ✅ |
| 2.2 | Local cohort enrichment (gnomAD MENA, KSA public) | 4h | High | ❌ |
| 2.3 | Variant prioritisation rules (ACMG-like filters) | 4h | High | ❌ |
| 2.4 | Clinical report generator (HTML, template) | 1 day | High | ✅ |
| 2.5 | IGV session generator + `localhost:60151` URLs | 2h | Medium | ✅ |
| 2.6 | Pangenome end-to-end integration test | 2h | High | ❌ |

---

## Phase 3 — Real data and benchmark (2-3 weeks)

Move from synthetic to real data, produce quantitative
evidence that the approach is better than linear.

| # | Item | Effort | Impact | Status |
|---|---|---|---|---|
| 3.1 | Real trio (KSA001 or public) end-to-end | 4h | High | ❌ |
| 3.2 | Benchmark: linear vs pangenome (F1-score) | 1 day | High | ❌ |
| 3.3 | Synthetic ground truth with Truvari | 1 day | High | ❌ |
| 3.4 | Pitch document for KFSH&RC / KAUST | 2h | High | ❌ |

---

## Phase 4 — Research and publication (4-6 weeks)

| # | Item | Effort | Impact | Status |
|---|---|---|---|---|
| 4.1 | HPO extraction without GPU (small model) | 2 days | High | ❌ |
| 4.2 | NL → CaseManifest agent (LLM parses clinical text) | 2 days | High | ❌ |
| 4.3 | Mini-tutorial / sandbox framework (parallel project) | 3 days | Medium | ❌ |
| 4.4 | Paper draft (Bioinformatics / GigaScience) | 3 days | High | ❌ |

---

## Fase 2.5 — Integración Sandbox (completada)

PublicGenomicAgent consume catálogos de casos sintéticos
producidos por el proyecto paralelo **Sandbox**.

### Contrato

- Formato: `manifest.yaml` + BAMs del trío + referencia
- Identificación por **SPDI canónico**
- Nombrado de directorio: SPDI con `:` → `_`
- Sin acoplamiento con el generador: PGA solo lee

### Implementado

| # | Item | Estado |
|---|---|---|
| 2.5.1 | Módulo `agent/sandbox.py` (SandboxCatalog + SandboxCase) | ✅ |
| 2.5.2 | Conversión a `CaseManifest` (pedigrí, ROIs, HPO) | ✅ |
| 2.5.3 | CLI `pga sandbox list` | ✅ |
| 2.5.4 | CLI `pga sandbox describe` | ✅ |
| 2.5.5 | CLI `pga sandbox run` (+ `--clinvar`, `--planner`) | ✅ |
| 2.5.6 | Tests unitarios (16) | ✅ |
| 2.5.7 | Test de integración end-to-end | ✅ |

### Verificado

```
Caso CYP2D6 (Sandbox)
  → pga sandbox run --reference test_ref.fa
  → 12 iteraciones OK
  → 1 variante recesiva detectada:
     chr2:42128941 A>T
     father 0/1, mother 0/1, proband 1/1
  → auto_rec_hom, score 3.0
```

### Pendiente (proyecto Sandbox)

- `ref.fa` en el contrato (autocontención del caso).
- Bug: `ground_truth.chrom` de CYP2D6 (`chr2` → `chr22`).
- Simular delecciones reales (CIGAR `D`) — NPHP1 no produce
  variantes detectables hoy.

### Pendiente (PGA)

- `pga sandbox run --clinvar` end-to-end (13 pasos).
- Tests de integración con `--planner llm`.

## Known issues / fixes pending

- **`LLMPlanner` puro (LLM-first)**: `_resolve_paths` y orden
  correcto del chequeo de duplicados. No usado en pipeline
  por defecto (usamos `RuleFirstLLMPlanner`). Fix estructural
  documentado en el commit `a224c97`.
- **`make_bam_with_variant.py`**: reporta `Lecturas con variante
  (CIGAR != full match): 0` aunque el pileup sí ve la variante.
  Fix cosmético del contador.
- **`docs/known_issues.md`**: confirmar explícitamente el fix
  BAQ ya aplicado en `call_variants.py` y `joint_call`.
- **Coverage**: no hay `pytest --cov` configurado. Medir y
  decidir umbral mínimo (objetivo: 80% en `src/`).
- **Lint / type check**: Ruff y mypy no están configurados.
  Añadir como `dev` deps de `pga-core`.

## Design decisions worth remembering

- **Reglas primero, LLM después**. `RuleFirstLLMPlanner`
  consulta al `RuleBasedPlanner` primero; solo llama al LLM
  si las reglas no saben qué hacer. Los modelos 3B no siguen
  instrucciones complejas como "no repitas tools".
- **Validación doble**. El LLM propone, Pydantic valida.
  Fallback determinista cuando algo falla.
- **Privacidad por arquitectura**. Basenames en el prompt,
  Ollama local por defecto, `include_paths=False`, claves de
  `state.vcfs` y `state.artifacts` como etiquetas lógicas.
- **Virtualización invisible**. El código nunca sabe que existe
  micromamba. Los entornos se resuelven vía `ToolRuntime`.
- **Trazabilidad completa**. Cada tool call queda registrada
  en `state.tool_calls` con input validado, output serializado,
  duración y error.

## How to resume work

```bash
cd ~/PublicGenomicAgent
git log --oneline -5
git status
pytest tests/ -q                       # 98 tests

# Sanity check
pga llm ping                           # local Ollama
pga env verify                         # 5 envs

# End-to-end demo
python3 scripts/make_demo_trio.py
pga run-trio \\
    --father results/demo_trio/father.bam \\
    --mother results/demo_trio/mother.bam \\
    --proband results/demo_trio/proband.bam \\
    --reference results/demo_trio/ref.fa \\
    --region chr1:1000-1200 \\
    --case-id DEMO_TRIO --label DEMO1
```

### Benchmark: 2-variant duplication (blocked)

The synthetic SD benchmark works with 1 variant (chr1:1050 T>G,
pangenome recovers it, linear misses it). Extending to 2 variants
is blocked: vg giraffe only aligns 41/100 reads to the graph,
regardless of paralog SNPs added to copy 2.

Hypothesis: the graph with two 200 bp copies + alt paths is too
large for giraffe with default parameters. Next steps:
  - Try vg giraffe with --num-em-batch or lower --min-identity.
  - Try with copies of 100 bp (as in the original benchmark)
    but variants separated by 30 bp instead of 100 bp.
  - Or accept the 1-variant benchmark as sufficient.

The 1-variant benchmark is already documented in
docs/benchmark_duplications.md and supports the pitch.
## [Día 2] Tool: `detect_homozygous_regions` (ROH restringido a panel)

**Motivo**: el diferenciador declarado del proyecto es soberanía + casos
consanguíneos. Un detector de regiones de homocigosidad aparente dentro
del panel convierte ese diferenciador en algo cuantificable y auditable.

**No es un ROH caller genómico.** Los paneles de genes no tienen densidad
de SNPs suficiente para estimar F_ROH global ni mapear ROH largos (>8 Mb).
La tool detecta "regiones de homocigosidad aparente" dentro de la zona
panelada, como heurística de priorización de variantes recesivas
homocigotas.

**Implementación prevista**:
- `src/publicgenomicagent/tools/roh.py` con Pydantic I/O.
- Wrapper de `bcftools roh` (no PLINK: más falsos positivos en panel).
- Parámetros conservadores: 2–3 het permitidos, longitud mínima 1–2 Mb.
- Cruce con BED del panel para anotar genes solapados.
- Campo `caveat` **obligatorio** en el output: "estimación restringida
  al panel, no interpretable como F_ROH genómico".
- Evento UI `roh_summary` + render propio (distinto de tool_result).
- Test unitario con VCF sintético (subprocess mockeado) + test de
  integración marcado `slow`.
- Escenario documentado en `docs/use_cases.md`: "consanguinidad sospechada".

**Pregunta abierta antes de implementar**: ¿el VCF de entrada trae `PL`
en FORMAT o solo `GT`? Determina si `bcftools roh` corre con likelihoods
o solo con genotipos llamados. Verificar con:
`bcftools query -f '%FORMAT\n' <vcf> | head`

**Bloqueante**: no empezar hasta que Día 1 esté commiteado y `pga chat`
arranque con Ollama.
## [Día 3+] Sembrar variantes reales en el caso demo

El caso demo (`cases/demo_trio/`) genera BAMs con reads aleatorios,
así que `call_variants` produce ~20 000 variantes ruidosas y
`mendelian_filter` las filtra todas (0 auto_rec_hom, 0 auto_dom, etc.).

Para que el caso sirva también para demos clínicas (pitch, Alkuraya),
`scripts/make_demo_case.py` debería plantar variantes concretas:
  - Una variante auto-recesiva homocigota en el proband.
  - Los padres HET en la misma posición.
  - Idealmente dentro de la ROI (chr2:40000-60000).

Eso permite verificar que `mendelian_filter` produce al menos 1
`auto_rec_hom_vcf` con contenido, y `prioritize_variants` puntúa algo.

Coste estimado: ~1 hora. Toca `make_bam()` en el script.
## [Día 3+] Actualizar actions de CI cuando migren a Node 24

GitHub ha deprecado Node 20. Los jobs de CI emiten warnings:
- actions/checkout@v4
- actions/setup-python@v5
- mamba-org/setup-micromamba@v1

Cuando salgan las versiones con Node 24, actualizar los `uses:` del
workflow. No urgente: hoy solo son warnings cosméticos.

También: `ubuntu-latest` migra a Ubuntu 26 en octubre de 2026.
Revisar que los entornos micromamba sigan resolviendo sin cambios.
## [Día 5+] Diagnosticar JSON truncado en nl_to_case_manifest

La tool `nl_to_case_manifest` falla con modelos pequeños (qwen2.5:3b
y 7b) porque el LLM devuelve un JSON truncado a mitad. Síntomas:

- El output empieza correctamente con '{' y los primeros campos.
- Se corta en medio de un string o de un objeto.
- `_extract_json_block` no encuentra un bloque balanceado.

Hipótesis a verificar:

1. **num_predict insuficiente**: subir a 3000 y probar de nuevo.
2. **El LLM se queda sin contexto**: `num_ctx` es 8192, el prompt
   ~1500 tokens, el output ~800-1000 tokens. Debería caber, pero
   Ollama puede reservar memoria de otra forma.
3. **El LLM añade texto antes/después del JSON** que rompe el parser:
   capturar la respuesta cruda y ver.
4. **El LLM emite `null` como string** en algunos campos (`affected: null`
   en lugar de `false`), lo que Pydantic rechazaría — pero el error
   actual es de parseo, no de validación.

Plan de diagnóstico:

- Añadir un flag `--debug` a `pga case parse` que imprima la respuesta
  cruda del LLM antes de parsear.
- O añadir un `logging.debug` en `nl_to_case_manifest` que guarde la
  respuesta a un fichero `raw_llm_response.txt` en el output_dir.

Coste estimado: 1-2 horas.

Prioridad: alta. Bloquea la utilidad real de la tool.

Nota sobre hardware:
- Con 15 GB RAM total, `qwen2.5:7b` + num_ctx:8192 no cabe y hace swap.
- `qwen2.5:3b` sí cabe, pero tarda ~2 min por manifest en CPU.
- Para el pitch: el pipeline escala el modelo según los recursos
  disponibles (3b en portátil, 7b+ en servidor con GPU)."""

## [Pendiente] Pulir bloque NL → CaseManifest

El pipeline `pga case parse` funciona end-to-end con qwen2.5:3b sobre
relatos en inglés. Quedan tres tareas de pulido:

1. **Falso positivo en warning de samples**: father/mother aparecen
   en el texto como "parents" / "padres" pero `_verify_samples_in_text`
   no los reconoce (SYNONYMS no cubre el plural). Añadir:
   - `father`: `["parents", "padres", "progenitors"]`
   - `mother`: `["parents", "padres", "progenitors"]`
   - `proband`: `["son", "daughter", "child", "male patient", "female patient"]`

2. **Probar con relato español**: solo probamos el inglés en la
   versión final. El prompt es bilingüe; falta verificar.

3. **Tests unitarios con mock del LLM**: no hay tests específicos
   para `nl_to_case_manifest` con el LLM mockeado. Hay tests de
   `_sanitize_lists`, `_validate_and_clean`, etc., pero no del
   flujo end-to-end. Añadir ~8-10 tests con `MagicMock` del cliente.

4. **Documentar hardware en README**:
   - Con 15 GB RAM: qwen2.5:3b funciona (~1-2 min/manifest).
   - Con 15 GB RAM: qwen2.5:7b hace swap, inviable.
   - Con >32 GB RAM o GPU: qwen2.5:7b va cómodo.
   - El pipeline escala el modelo según los recursos disponibles.

Coste estimado: ~1 hora.

Prioridad: media. Bloquea parcialmente la calidad del NL, pero no
el resto del pipeline.

## [Pendiente] Email a Alkuraya

Documento de presentación para Fowzan Alkuraya (KFSHRC) proponiendo
colaboración en desarrollo conjunto. Munición técnica disponible
hoy:

- Caso demo con variante recesiva sembrada (padres HET, proband HOM).
- GIAB HG002 real procesado end-to-end (3 tools, 1765 variantes).
- `detect_homozygous_regions` (ROH restringido a VCF).
- `nl_to_case_manifest` (relato clínico → manifest estructurado).
- `pga audit` para verificación de integridad de sesiones.
- CI verde con 185 tests + smoke E2E.

Pendiente:
- Redactar el cuerpo del email (2 páginas).
- Preparar el apéndice técnico (1-2 páginas).
- Adjuntos: session.json del caso GIAB, session.json del caso demo,
  fragmento de `pga audit --verbose`.
- Opcional: vídeo de 2 min con asciinema del pipeline corriendo.

Coste estimado: 2-3 horas.

Prioridad: alta. Es la visibilidad externa del proyecto.
## [Día 6+] Manifest-first enforcement: no ejecutar sin manifest completo y correcto

**Principio de diseño**: un `CaseManifest` completo y correcto es la
única puerta de entrada al pipeline. No hay atajos. Si falta algo,
el sistema lo pide activamente antes de ejecutar cualquier tool.

**Por qué**: hoy el sistema falla a mitad del pipeline si falta un
campo (`CaseManifest.validate_minimum()` o el `ValueError` de una
tool). Eso obliga al usuario a conocer la gramática del workflow y
a reiniciar. Un clínico no debería aprender esa gramática: debería
describir el caso y el sistema debería pedir lo que falte.

### Flujo objetivo

1. Usuario: `pga chat`
2. Sistema: "¿Puedes describir el caso?"
3. Usuario: relato libre (EN/ES)
4. Sistema (LLM): extrae manifest preliminar.
5. Sistema (validación): verifica cada campo contra el texto original.
6. Sistema: "El manifest está incompleto. ¿Qué genes o regiones
   quieres analizar?"
7. Usuario: "BRCA1, BRCA2"
8. Sistema (validación): verifica que "BRCA1" y "BRCA2" aparecen en
   la conversación.
9. Sistema: "Voy a analizar BRCA1 y BRCA2 en un trío (proband
   afectado, padres sanos, no consanguíneos). ¿Correcto? [y/N]"
10. Usuario: "y"
11. Sistema: ejecuta el pipeline. Registra todo en session.json.

En cada paso, el sistema no avanza sin validar y sin confirmar.

### Implementación

**Bloque 1 — Validación (medio día)**
- `src/publicgenomicagent/tools/case_manifest.py`:
  - `required_fields_for_pipeline(manifest) -> dict[str, list[str]]`
    devuelve qué campos faltan para cada tool del pipeline por
    defecto.
  - `validate_manifest_against_source(manifest, text) -> list[str]`
    verifica que cada campo del manifest está respaldado por el
    texto original (samples, genes, consanguinidad, ROIs).
- Tests unitarios de ambas funciones.

**Bloque 2 — Diálogo (1 día)**
- `src/publicgenomicagent/ui/chat.py`:
  - `ensure_manifest_ready(manifest, ui) -> CaseManifest`: si faltan
    campos, lanza preguntas al usuario y actualiza el manifest.
  - Extender el REPL para hacer preguntas multi-turno.
- `src/publicgenomicagent/agent/loop.py`:
  - Antes de `_run_one`, verificar que el manifest está listo para
    la tool que se va a ejecutar. Si no, no ejecutar y devolver un
    evento `needs_user_input`.

**Bloque 3 — Integración y tests (medio día)**
- Integrar todo en `pga chat` y `pga case parse`.
- Tests de comportamiento:
  - Manifest completo → pipeline arranca sin preguntas.
  - Manifest sin ROIs → pregunta por genes/regiones.
  - Manifest sin pedigree → pregunta por trío/familia.
  - LLM propone un gen que no está en el texto → rechazado.
  - Usuario no confirma → pipeline no arranca.

**Total estimado**: 2-3 días.

### Consecuencias en el pitch

El pitch debe decir:

> "PublicGenomicAgent no ejecuta ningún análisis sin un CaseManifest
> completo y correcto. Si el usuario no ha proporcionado toda la
> información, el sistema la pide activamente, valida que cada campo
> esté respaldado por evidencia del input, y solo ejecuta el pipeline
> cuando el usuario confirma explícitamente. El LLM propone, el
> sistema valida, el usuario confirma. No hay atajos."

### Lo que NO hace

- No ejecuta "un poco" y completa después. O hay manifest, o no hay
  análisis.
- No acepta campos inventados por el LLM sin evidencia en el texto.
- No avanza sin confirmación explícita del usuario.

### Prioridad

Alta. Bloquea el uso clínico real. Pero no urgente para el email
a Alkuraya: se puede mencionar como principio de diseño aunque la
implementación completa llegue después.
