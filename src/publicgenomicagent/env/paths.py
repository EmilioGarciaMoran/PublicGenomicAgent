from __future__ import annotations

import os
from pathlib import Path


PGA_ROOT = Path(os.environ.get("PGA_ROOT", Path.home() / ".pga"))
ENVS_DIR = PGA_ROOT / "envs"
CACHE_DIR = PGA_ROOT / "cache"
LOGS_DIR = PGA_ROOT / "logs"
LOCKS_DIR = PGA_ROOT / "locks"


def ensure_pga_root() -> None:
    for d in (PGA_ROOT, ENVS_DIR, CACHE_DIR, LOGS_DIR, LOCKS_DIR):
        d.mkdir(parents=True, exist_ok=True)


def repo_root() -> Path:
    """Raíz del repo. Asume src/publicgenomicagent/env/paths.py."""
    return Path(__file__).resolve().parents[3]


def registry_path() -> Path:
    """Localiza envs/registry.yaml desde cwd o desde el paquete."""
    candidates = [
        Path.cwd() / "envs" / "registry.yaml",
        repo_root() / "envs" / "registry.yaml",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("envs/registry.yaml no encontrado")
