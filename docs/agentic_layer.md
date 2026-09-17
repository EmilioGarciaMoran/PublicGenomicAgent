
# Capa agéntica (recordatorio de diseño)

Documento breve para fijar decisiones antes de construir el agente LLM.
No es diseño exhaustivo, es una restricción y un mapa.

## Regla de diseño base

Toda tool debe ser:

- Función pura, sin estado global.
- Entrada y salida tipadas con Pydantic.
- Argumentos serializables a JSON.
- Invocable desde Python sin pasar por CLI.
- Auditable: cada ejecución registra tool, args y resultado.

Sin esto, el agente no podrá invocar tools de forma fiable.

## Niveles de auditoría del agente

1. Generación de prompts
   - Baterías de prompts: clínicos, ambiguos, informales, coordenadas directas.
   - Verificar que el agente infiere correctamente la tool a invocar.

2. Benchmark de selección de tools
   - Matriz input vs tools esperadas.
   - Validación de argumentos.
   - Herramientas prohibidas por caso.

3. Inyección de fallos
   - BAM corrupto, índice ausente.
   - VCF vacío.
   - Coordenadas fuera de rango.
   - Sample name inexistente.
   - Entorno micromamba no creado.
   - Verificar mensajes útiles, sin bucles ni tracebacks.

## Estructura prevista

tests/evals/
  TC_NEPHRO_01.json
  TC_BRUGADA_01.json
  ...

Formato mínimo de cada eval:

  {
    "id": "TC_NEPHRO_01",
    "input_prompt": "...",
    "expected_tool_calls": [
      {"tool": "fetch_roi", "args": {...}},
      {"tool": "compare_vcfs", "args": {...}}
    ],
    "forbidden_tools": ["..."],
    "mocks": { "call_variants": {"return": {...}} }
  }

Sin código: los evals son datos. Se añaden casos sin tocar Python.

## Mocking

Las tools del sistema deben poder sustituirse por mocks sin
reescribir el agente. Formalizar en tools/registry.py:
registro de implementación real y de mock, intercambiables.

## Estado

Pendiente. Se implementa tras local_pangenome y mendelian.
Este documento solo fija la restricción de diseño para no
construir tools incompatibles con el agente.
