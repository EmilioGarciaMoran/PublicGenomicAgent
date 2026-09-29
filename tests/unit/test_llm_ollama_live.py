"""Test end-to-end contra Ollama real.

Se salta automáticamente si el demonio Ollama no responde en
http://localhost:11434. Pensado para correr en el MacAir del
usuario, no en CI.
"""
from __future__ import annotations

import json
import urllib.request

import pytest

from publicgenomicagent.agent.llm import LLMError, OllamaClient
from publicgenomicagent.agent.llm_factory import build_llm_client
from publicgenomicagent.agent.config import LLMConfig


OLLAMA_TAGS = "http://localhost:11434/api/tags"


def _ollama_available() -> bool:
    try:
        with urllib.request.urlopen(OLLAMA_TAGS, timeout=2) as r:
            data = json.loads(r.read().decode("utf-8"))
        return len(data.get("models", [])) > 0
    except Exception:  # noqa: BLE001
        return False


def _first_model() -> str | None:
    try:
        with urllib.request.urlopen(OLLAMA_TAGS, timeout=2) as r:
            data = json.loads(r.read().decode("utf-8"))
        models = data.get("models", [])
        return models[0]["name"] if models else None
    except Exception:  # noqa: BLE001
        return None


pytestmark = pytest.mark.skipif(
    not _ollama_available(),
    reason="Ollama no disponible en localhost:11434",
)


def test_ollama_live_completes_a_trivial_prompt():
    model = _first_model()
    assert model is not None

    client = OllamaClient(
        endpoint="http://localhost:11434/api/generate",
        model=model,
        timeout=120,
    )
    out = client.complete(
        system='Responde SOLO con el JSON {"ok": true}.',
        user="Devuelve el JSON pedido.",
    )
    assert isinstance(out, str)
    assert len(out) > 0


def test_factory_builds_ollama_client_from_config():
    client = build_llm_client(
        LLMConfig(
            provider="ollama",
            endpoint="http://localhost:11434/api/generate",
            model=_first_model() or "qwen2.5:7b",
            timeout_seconds=120,
        )
    )
    assert isinstance(client, OllamaClient)
    out = client.complete("sys", "ping")
    assert isinstance(out, str)
