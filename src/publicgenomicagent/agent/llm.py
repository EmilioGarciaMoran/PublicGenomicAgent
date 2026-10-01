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


# ---------------------------------------------------------------------------
# Clientes HTTP (sin SDKs, solo urllib de la stdlib)
# ---------------------------------------------------------------------------

import json as _json
import urllib.error as _urlerror
import urllib.request as _urlreq


class HttpLLMClient:
    """Cliente LLM vía HTTP genérico.

    El contrato es mínimo: `complete(system, user) -> str`. Cada
    subclase concreta (Ollama, llama.cpp, ...) solo tiene que
    sobrescribir `_build_payload` y `_extract_text`.

    Usamos `urllib` de la stdlib a propósito: cero dependencias
    nuevas para un cliente que debe poder correr en un MacAir sin
    red, sin `pip install` adicionales.
    """

    def __init__(
        self,
        endpoint: str,
        model: str,
        timeout: float = 120.0,
        temperature: float = 0.0,
    ):
        self.endpoint = endpoint
        self.model = model
        self.timeout = timeout
        self.temperature = temperature

    # --- API pública -----------------------------------------------------

    def complete(self, system: str, user: str) -> str:
        payload = self._build_payload(system, user)
        body = _json.dumps(payload).encode("utf-8")
        req = _urlreq.Request(
            self.endpoint,
            data=body,
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        try:
            with _urlreq.urlopen(req, timeout=self.timeout) as resp:
                raw = resp.read().decode("utf-8")
        except _urlerror.URLError as e:
            raise LLMError(f"error de red con {self.endpoint}: {e}") from e
        except TimeoutError as e:
            raise LLMError(f"timeout tras {self.timeout}s") from e

        try:
            data = _json.loads(raw)
        except _json.JSONDecodeError as e:
            raise LLMError(f"respuesta no JSON: {raw[:200]}") from e

        return self._extract_text(data)

    # --- Hooks para subclases -------------------------------------------

    def _build_payload(self, system: str, user: str) -> dict:
        raise NotImplementedError

    def _extract_text(self, data: dict) -> str:
        raise NotImplementedError


class OllamaClient(HttpLLMClient):
    """Cliente para Ollama (`/api/generate`)."""

    def __init__(
        self,
        endpoint: str,
        model: str,
        timeout: float = 120.0,
        temperature: float = 0.0,
        num_ctx: int = 2048,
    ):
        super().__init__(
            endpoint=endpoint,
            model=model,
            timeout=timeout,
            temperature=temperature,
        )
        # Contexto efectivo de Ollama. El prompt real de PublicGenomicAgent
        # ronda los 800-1000 tokens; 2048 deja margen y evita que Ollama
        # reserve 32K por defecto, lo cual penaliza la latencia y la RAM.
        self.num_ctx = num_ctx

    def _build_payload(self, system: str, user: str) -> dict:
        return {
            "model": self.model,
            "system": system,
            "prompt": user,
            "stream": False,
            # format=json fuerza a Ollama a restringir la salida a JSON
            # válido. Sin esto, Qwen 7B suele devolver markdown o texto
            # libre, lo cual obliga a fallback al RuleBasedPlanner.
            "format": "json",
            "options": {
                "temperature": self.temperature,
                "num_ctx": self.num_ctx,
                "num_predict": 128,
            },
        }

    def _extract_text(self, data: dict) -> str:
        if "response" not in data:
            raise LLMError(
                f"respuesta de Ollama sin campo 'response': {list(data)[:5]}"
            )
        return data["response"]


class LlamaCppClient(HttpLLMClient):
    """Cliente para `llama-server` (`/completion`)."""

    def _build_payload(self, system: str, user: str) -> dict:
        # llama.cpp no tiene campo system nativo: lo concatenamos.
        prompt = f"{system}\n\n{user}" if system else user
        return {
            "prompt": prompt,
            "temperature": self.temperature,
            "n_predict": 1024,
            "stream": False,
        }

    def _extract_text(self, data: dict) -> str:
        if "content" not in data:
            raise LLMError(
                f"respuesta de llama.cpp sin campo 'content': {list(data)[:5]}"
            )
        return data["content"]


class NullLLMClient:
    """Cliente que siempre falla. Útil para forzar fallback en tests."""

    def complete(self, system: str, user: str) -> str:
        raise LLMError("NullLLMClient: no hay backend LLM configurado")
