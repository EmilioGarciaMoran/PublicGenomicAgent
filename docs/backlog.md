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
