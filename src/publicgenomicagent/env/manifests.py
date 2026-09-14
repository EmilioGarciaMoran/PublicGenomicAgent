from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .paths import LOCKS_DIR


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def lock_path(env_name: str) -> Path:
    return LOCKS_DIR / f"{env_name}.lock.json"


def read_lock(env_name: str) -> dict | None:
    p = lock_path(env_name)
    if not p.exists():
        return None
    return json.loads(p.read_text())


def write_lock(env_name: str, data: dict) -> None:
    p = lock_path(env_name)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2, sort_keys=True))
