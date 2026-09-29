"""Clientes LLM inyectables.

Por defecto, el proyecto NO habla con ningún proveedor externo. El
`FakeLLMClient` permite testear el `LLMPlanner` sin red, sin API key
y sin dependencias nuevas.

Para producción, se implementa un `LLMClient` que envuelva el SDK
del proveedor elegido (OpenRouter, Anthropic, OpenAI, Mistral, o
un servidor local con Ollama/vLLM). El contrato es deliberadamente
mínimo: `complete(system, user) -> str`. No hay streaming, tool
calling nativo, ni estado entre llamadas. Esto es intencional:

  - El catálogo de tools es cerrado (TOOL_REGISTRY). El LLM no
    decide llamar a URLs arbitrarias.
  - El planner es quien valida, no el cliente.
  - Un contrato mínimo es fácil de mockear y fácil de auditar.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Protocol


class LLMClient(Protocol):
    """Contrato mínimo de un cliente LLM."""

    def complete(self, system: str, user: str) -> str:
        """Devuelve la respuesta cruda del modelo como texto."""
        ...


@dataclass
class LLMResponse:
    """Respuesta tipada del planner (post-parseo del JSON del LLM)."""

    tool_name: str | None
    args: dict[str, Any]
    rationale: str
    stop: bool = False
    raw: str = ""


class LLMError(RuntimeError):
    """Error controlado del cliente LLM."""


class FakeLLMClient:
    """Cliente LLM determinista para tests.

    Se le pasa una lista de respuestas (strings). Cada llamada a
    `complete` devuelve la siguiente. Si se agotan, lanza `LLMError`.
    """

    def __init__(self, responses: list[str]):
        self._responses = list(responses)
        self._index = 0
        self.calls: list[dict[str, str]] = []

    def complete(self, system: str, user: str) -> str:
        self.calls.append({"system": system, "user": user})
        if self._index >= len(self._responses):
            raise LLMError(
                f"FakeLLMClient agotado tras {self._index} llamadas"
            )
        resp = self._responses[self._index]
        self._index += 1
        return resp


@dataclass
class RecordingLLMClient:
    """Cliente que devuelve siempre la misma respuesta y graba los prompts.

    Útil para inspeccionar qué se le envió al LLM en cada test.
    """

    response: str
    calls: list[dict[str, str]] = field(default_factory=list)

    def complete(self, system: str, user: str) -> str:
        self.calls.append({"system": system, "user": user})
        return self.response


# ---------------------------------------------------------------------------
# Parseo de la respuesta JSON del LLM
# ---------------------------------------------------------------------------

class LLMResponseFormatError(ValueError):
    """La respuesta del LLM no es un JSON válido con la forma esperada."""


def parse_llm_response(raw: str) -> LLMResponse:
    """Parsea la respuesta del LLM como JSON estricto.

    Acepta dos formas:

        {"tool_name": str, "args": {...}, "rationale": str}
        {"stop": true, "rationale": str}

    No acepta markdown, bloques ```json, ni texto adicional. El
    SYSTEM_PROMPT ya pide JSON puro; si el LLM no lo cumple, es
    responsabilidad del planner hacer fallback.
    """
    text = raw.strip()
    if not text:
        raise LLMResponseFormatError("respuesta vacía")

    # Tolerancia mínima: si viene envuelto en ```json ... ```, lo pelamos.
    # Es común que los LLM añadan fences aunque se les pida JSON puro.
    if text.startswith("```"):
        lines = text.splitlines()
        # Quitar primera línea (``` o ```json) y última (```)
        if lines and lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].strip().startswith("```"):
            lines = lines[:-1]
        text = "\n".join(lines).strip()

    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        raise LLMResponseFormatError(f"JSON inválido: {e}") from e

    if not isinstance(data, dict):
        raise LLMResponseFormatError(
            f"se esperaba un objeto JSON, se obtuvo {type(data).__name__}"
        )

    if data.get("stop") is True:
        return LLMResponse(
            tool_name=None,
            args={},
            rationale=str(data.get("rationale", "")),
            stop=True,
            raw=raw,
        )

    tool_name = data.get("tool_name")
    if not isinstance(tool_name, str) or not tool_name:
        raise LLMResponseFormatError(
            "falta 'tool_name' (string no vacío) o 'stop'=true"
        )

    args = data.get("args", {})
    if not isinstance(args, dict):
        raise LLMResponseFormatError("'args' debe ser un objeto JSON")

    rationale = str(data.get("rationale", ""))

    return LLMResponse(
        tool_name=tool_name,
        args=args,
        rationale=rationale,
        stop=False,
        raw=raw,
    )
