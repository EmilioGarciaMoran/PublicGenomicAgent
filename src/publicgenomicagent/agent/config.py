"""Configuración del agente: LLM, privacidad, rutas de caché.

Lee `~/.pga/config.yaml`. Si no existe, escribe uno con defaults
y lo devuelve. Nunca escribe secretos: las API keys (si algún día
se usan) deben ir por variables de entorno, no por este fichero.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from ..env.paths import PGA_ROOT


CONFIG_PATH = PGA_ROOT / "config.yaml"


DEFAULT_CONFIG: dict = {
    "llm": {
        "provider": "ollama",  # ollama | llamacpp | none
        "endpoint": "http://localhost:11434/api/generate",
        "model": "qwen2.5:7b",
        "timeout_seconds": 120,
        "temperature": 0.0,
        "include_paths": False,  # privacidad por defecto
        "num_ctx": 2048,
    }
}


@dataclass
class LLMConfig:
    provider: str = "ollama"
    endpoint: str = "http://localhost:11434/api/generate"
    model: str = "qwen2.5:7b"
    timeout_seconds: float = 120.0
    temperature: float = 0.0
    include_paths: bool = False
    num_ctx: int = 2048


@dataclass
class AgentConfig:
    llm: LLMConfig = field(default_factory=LLMConfig)
    raw: dict = field(default_factory=dict)


def _write_default(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(DEFAULT_CONFIG, sort_keys=False))


def load_config(path: Path | None = None) -> AgentConfig:
    """Carga la config del usuario. Si no existe, la crea con defaults."""
    p = path or CONFIG_PATH
    if not p.exists():
        _write_default(p)

    raw = yaml.safe_load(p.read_text()) or {}
    llm_raw = raw.get("llm", {})

    llm = LLMConfig(
        provider=llm_raw.get("provider", "ollama"),
        endpoint=llm_raw.get("endpoint", "http://localhost:11434/api/generate"),
        model=llm_raw.get("model", "qwen2.5:7b"),
        timeout_seconds=float(llm_raw.get("timeout_seconds", 120)),
        temperature=float(llm_raw.get("temperature", 0.0)),
        include_paths=bool(llm_raw.get("include_paths", False)),
        num_ctx=int(llm_raw.get("num_ctx", 2048)),
    )

    # Variables de entorno pueden sobreescribir (útil para tests y CI).
    if v := os.environ.get("PGA_LLM_ENDPOINT"):
        llm.endpoint = v
    if v := os.environ.get("PGA_LLM_MODEL"):
        llm.model = v

    return AgentConfig(llm=llm, raw=raw)
