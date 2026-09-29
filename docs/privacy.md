# Política de privacidad de PublicGenomicAgent

Este documento formaliza las garantías de privacidad del sistema.
No es una declaración de intenciones: cada afirmación está respaldada
por código y por tests unitarios reproducibles.

## Principios

1. **El paciente no sale del proceso.** Ningún dato identificativo
   del paciente abandona la máquina local salvo que el operador
   active explícitamente un cliente LLM remoto. Por defecto, el
   sistema usa un LLM local (Ollama en `localhost:11434`).
2. **Las rutas absolutas no se envían al LLM.** El prompt que recibe
   el planificador contiene solo basenames (`proband.bam`), no rutas
   completas (`/home/hospital/paciente_1234.exoma.bam`). Esto es
   deliberado: los directorios clínicos suelen incluir identificadores.
3. **Sin red en el camino crítico.** El análisis genómico
   (samtools, bcftools, vg, plink) se ejecuta contra binarios locales
   instalados en entornos micromamba aislados. No hay llamadas a
   servicios externos en el flujo de análisis.
4. **Auditoría completa.** Cada decisión del agente queda registrada
   en `AgentState.tool_calls` con input validado, output serializado,
   duración y error. El estado completo se persiste en `session.json`
   para reproducibilidad y para auditoría posterior.
5. **Sin secretos en ficheros de configuración.** `~/.pga/config.yaml`
   no contiene API keys. Si en el futuro se añade un proveedor remoto,
   las credenciales irán por variables de entorno, nunca por el YAML.

## Qué se envía al LLM

El planificador (`LLMPlanner`) construye el prompt con
`render_user_prompt(state, include_paths=False)` por defecto.

**Se envía:**

- `case_id`, `source`, `proband`, `consanguinity`
- Pedigrí: samples, sexo, afectación, relaciones padre/madre
- ROIs candidatos: `chrom:start-end` y etiqueta (nombre del gen)
- Términos HPO por sample, si los hay
- **Basenames** de los recursos: `proband.bam`, `hg38.fa`, `cohort.vcf.gz`
- Historial resumido de tool calls: `fetch_roi: ok`, `qc_bam: error=...`
- Catálogo de tools: nombre, descripción, campos de input

**NO se envía:**

- Rutas absolutas (`/home/...`, `/mnt/...`, `/data/...`)
- Contenido de los BAM, VCF, FASTQ
- Secuencias genómicas
- Identificadores del paciente que no sean los `sample` del pedigrí
- Cualquier metadato del sistema operativo o del usuario

El operador puede activar `--include-paths` o `include_paths: true`
en la config **solo si el LLM es estrictamente local y el operador
asume la responsabilidad**. Esto está documentado como excepción,
no como default.

## Arquitectura de privacidad por capas

```
[ Operador ]
     |
     v
[ CLI / Chat UI ]
     |
     v
[ AgentLoop + LLMPlanner ]  -- envía SOLO basenames + metadatos -->  [ LLM local ]
     |                                                                  (Ollama / llama.cpp)
     v
[ SessionContext.call_tool ]
     |
     v
[ ToolRuntime ]  -- resuelve binarios en entornos micromamba locales -->  [ samtools, bcftools, vg, plink ]
     |
     v
[ Sistema de ficheros local ]  -- BAMs, VCFs, referencias del usuario -->  [ nunca salen de aquí ]
```

El único punto de salida de información desde el proceso es el
`LLMClient`. Si el cliente es `NullLLMClient` (config `provider: none`),
el `LLMPlanner` cae al `RuleBasedPlanner` y **cero información sale
del proceso**.

## Garantías verificadas por tests

Cada una de estas afirmaciones está cubierta por tests unitarios:

| Garantía | Test |
|---|---|
| El prompt no contiene rutas absolutas | `test_llm_planner_does_not_send_paths_by_default` |
| El operador puede activarlas explícitamente | `test_llm_planner_sends_paths_when_explicitly_enabled` |
| Un cliente nulo fuerza fallback determinista | `test_null_client_always_raises` |
| El planificador valida antes de ejecutar | `test_llm_planner_falls_back_on_unknown_tool` |
| Args inválidos no llegan a la tool | `test_dispatch_invalid_input_raises_validation_error` |
| Cada tool call queda trazada | `test_record_appends_and_indexes` |

## Configuración recomendada para entornos clínicos

Para uso con datos de pacientes reales (hospital, KSA, UE), la
configuración debe ser:

```yaml
# ~/.pga/config.yaml
llm:
  provider: ollama
  endpoint: http://localhost:11434/api/generate
  model: qwen2.5:7b
  include_paths: false        # OBLIGATORIO en clínico
  num_ctx: 2048
```

Reglas adicionales:

- Ollama corriendo en la misma máquina que el análisis. Nunca
  `http://<host-remoto>:11434`.
- Sin conexión a internet durante el análisis, si es posible.
- `~/.pga/config.yaml` con permisos `600` (solo el usuario).
- BAMs, VCFs y referencias en disco cifrado (FileVault, LUKS).
- Los `session.json` generados contienen solo basenames y metadatos
  del caso, pero **deben tratarse como datos clínicos** y seguir la
  política de retención de la institución.

## Qué NO hace el sistema

Para evitar expectativas incorrectas:

- **No anonimiza los datos.** La privacidad viene de no enviarlos,
  no de transformarlos. Si el operador activa un proveedor remoto
  con `include_paths: true` y `session.json` con sample IDs reales,
  los datos salen tal cual.
- **No cifra los ficheros intermedios.** Los sub-BAMs, VCFs y
  grafos se escriben en `results/` sin cifrado. Depende del sistema
  de ficheros del operador.
- **No gestiona consentimiento informado ni cumplimiento normativo.**
  Eso es responsabilidad de la institución. El sistema solo garantiza
  que, por diseño, la ruta por defecto no expone datos fuera de la
  máquina.
- **No controla el LLM.** Si el operador configura un Ollama remoto
  o un proveedor en la nube, el sistema no lo bloquea. Confía en la
  configuración y en la auditoría posterior.

## Referencias

- `src/publicgenomicagent/agent/prompts.py` — `render_state()` y
  `render_tools_catalog()` con el comportamiento de `include_paths`.
- `src/publicgenomicagent/agent/llm.py` — `HttpLLMClient` (urllib
  stdlib, sin telemetría), `NullLLMClient`.
- `src/publicgenomicagent/agent/llm_factory.py` — resolución de
  `provider` a cliente.
- `tests/unit/test_planner.py` — tests de privacidad por defecto.
- `tests/unit/test_llm_clients.py` — tests contra `FakeHttpServer`,
  sin red externa.
- `docs/virtualization.md` — aislamiento por entornos micromamba.

