# NEXT_CHAT — Mensaje de apertura para el próximo chat

> **Instrucciones de uso**: copia TODO el bloque entre las líneas
> `=== INICIO MENSAJE ===` y `=== FIN MENSAJE ===` como primer
> mensaje del próximo chat con el asistente. No añadas nada más.

=== INICIO MENSAJE ===

Retomamos PublicGenomicAgent. Sprint 0, Día 1.

## Contexto del proyecto

**Objetivo**: agente genómico conversacional OFFLINE para tests genéticos,
ejecutable en un MacBook Air, multilingüe (EN/ES/AR), que produce un
informe auditable de asociación genotipo–fenotipo. Candidato a validación
clínica por Fowzan Alkuraya (KFSHRC, Arabia Saudí).

**Diferencial**: pangenoma local + soberanía de datos + auditabilidad
completa de procesos. No competimos en precisión bruta con Illumina/
DRAGEN: competimos en soberanía, transparencia y casos consanguíneos.

## Estado verificado (2026-10-08)

- Commit: `91b75b1` — "docs(readme): add Benchmarks section"
- Tests: **110 unit passed** (~20 s)
- LLM: **Ollama + Qwen2.5:3b ya integrado** (`agent/llm.py`, `agent/llm_factory.py`)
- Loop agéntico: `agent/loop.py` OK
- Planner: `agent/planner.py` OK
- Prompts: `agent/prompts.py` OK
- Informe HTML + IGV XML: OK (`results/DEMO_TRIO.report.html`, `.igv.xml`)
- Sandbox sintético + benchmarks: OK
- 13 tools en `TOOL_REGISTRY` con Pydantic I/O (`src/publicgenomicagent/tools/`)
- setup.sh idempotente pero **asume micromamba previo**
- CI: solo `ubuntu-latest`, integration con `continue-on-error: true`

## Gaps confirmados (4 pilares)

1. **INTERFAZ**: no hay UI chat; solo CLI determinista.
2. **TESTEO**: CI solo en ubuntu; integration tests no bloquean merge;
   no hay job de instalación limpia ni smoke end-to-end.
3. **DOCUMENTACIÓN**: install requiere micromamba previo; falta
   whitepaper para Alkuraya; README sin badges.
4. **AUDIT**: `session.json` existe pero sin trazas completas
   (hashes, versiones de tools, comandos reproducibles) ni comando
   `pga audit`.

## Decisión ya tomada (no deliberar más)

- **LLM**: Ollama + `qwen2.5:3b` (Q4). Ya está integrado.
  No cambiar. No evaluar alternativas en este chat.

## [X] del Día 1

**Construir la UI chat TUI sobre `agent/loop.py`, con entry point
`pga chat` y soporte de sesión persistente.**

Objetivo concreto del Día 1:
- `pga chat` arranca, muestra prompt, envía input al LLM local,
  recibe `tool_calls`, ejecuta las tools del registry, muestra
  respuesta final, permite varios turnos.
- Sesión persistente reutilizando `agent/state.py` + `session.json`.
- Sin i18n todavía (eso es Día 3).
- Sin informe HTML todavía (eso viene después).

## Lo que quiero del asistente en la primera respuesta

1. **Diseño de `src/publicgenomicagent/ui/`**: estructura de archivos
   y responsabilidades de cada uno.
2. **Decisión razonada** entre `Textual`, `prompt_toolkit` o solo `Rich`
   para esta TUI, con justificación corta (una o dos frases por opción).
3. **Esqueleto mínimo** de `pga chat` que arranca, habla con
   `agent/loop.py` y cierra limpiamente. Código real, no pseudocódigo.
4. **Plan de test mínimo** para ese Día 1 (qué probar, dónde, cómo).
5. **Commit sugerido** para cerrar Día 1.

No quiero roadmap a 3 semanas en la primera respuesta. Quiero
**Día 1 ejecutable**, con archivos concretos y comandos concretos.

## Documentos de referencia en el repo

- `docs/HANDOFF.md` — estado consolidado
- `docs/pitch.md` — pitch formal (3 pilares, benchmarks, roadmap)
- `docs/agentic_layer.md` — diseño de la capa agéntica
- `docs/quickstart.md` — walkthrough end-to-end
- `docs/use_cases.md` — 4 escenarios clínicos
- `docs/backlog.md` — backlog
- `docs/installation.md` — instalación (a reescribir)
- `/tmp/pga_resume_20261008.txt` — resume de la sesión anterior

## Primer paso concreto que te pido

Antes de escribir código, **pregúntame lo estrictamente necesario
para no equivocarte en el diseño de la UI**. Idealmente 3 preguntas o menos.
Si no necesitas preguntar nada, adelante y propón el diseño directamente.

Después, dame el esqueleto y arrancamos.

=== FIN MENSAJE ===

---

## Notas para el asistente (no copiar al chat)

1. **No reabrir la decisión del LLM.** Ya está tomada: Ollama + Qwen2.5:3b.
2. **No proponer Docker.** El objetivo es MacBook Air nativo, sin contenedores.
3. **No proponer servicios cloud.** El diferencial es la soberanía.
4. **No añadir features fuera del [X] del día.** El multi-idioma es Día 3.
5. **Respetar el stack existente**: Typer, Rich, Pydantic, pytest.
6. **Textual es la apuesta natural** para la TUI porque ya usas Rich.
   Confirmar antes de instalar dependencias nuevas.
7. **Cada respuesta, un único commit sugerido.**

## Checklist de cierre del Día 1

git add -A
git commit -m "feat(ui): add 'pga chat' TUI over agent/loop.py"
git push origin main

# Regenerar resume
{
  echo "=== GIT LOG ==="
  git log --oneline -10
  echo ""
  echo "=== TESTS ==="
  /home/egarmo/bin/micromamba run -p ~/.pga/envs/pga-core \
    pytest tests/unit/ -q --no-header 2>&1 | tail -5
} > /tmp/pga_resume_$(date +%Y%m%d).txt

# Actualizar docs/HANDOFF.md con la fecha del día
