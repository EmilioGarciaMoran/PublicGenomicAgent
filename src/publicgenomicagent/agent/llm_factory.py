"""Construcción del cliente LLM a partir de la configuración.

Reglas:

  - provider=ollama  -> OllamaClient
  - provider=llamacpp -> LlamaCppClient
  - provider=none    -> NullLLMClient (fuerza fallback al planner)
"""
from __future__ import annotations

from .config import LLMConfig
from .llm import (
    LLMClient,
    LlamaCppClient,
    NullLLMClient,
    OllamaClient,
)


def build_llm_client(cfg: LLMConfig) -> LLMClient:
    provider = (cfg.provider or "none").lower()
    if provider == "ollama":
        return OllamaClient(
            endpoint=cfg.endpoint,
            model=cfg.model,
            timeout=cfg.timeout_seconds,
            temperature=cfg.temperature,
            num_ctx=cfg.num_ctx,
        )
    if provider == "llamacpp":
        return LlamaCppClient(
            endpoint=cfg.endpoint,
            model=cfg.model,
            timeout=cfg.timeout_seconds,
            temperature=cfg.temperature,
        )
    return NullLLMClient()
