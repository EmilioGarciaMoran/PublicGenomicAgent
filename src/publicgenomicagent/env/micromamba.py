from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from .paths import ENVS_DIR, ensure_pga_root


def find_micromamba() -> str:
    candidates = [
        os.environ.get("MICROMAMBA_BIN"),
        str(Path.home() / "bin" / "micromamba"),
        shutil.which("micromamba"),
    ]
    for c in candidates:
        if c and Path(c).exists() and os.access(c, os.X_OK):
            return c
    raise RuntimeError(
        "micromamba no encontrado. Instálalo o exporta MICROMAMBA_BIN=/ruta/a/micromamba"
    )


def env_prefix(name: str) -> Path:
    return ENVS_DIR / name


def env_exists(name: str) -> bool:
    return (env_prefix(name) / "conda-meta").is_dir()


def bin_path(env_name: str, binary: str) -> Path:
    return env_prefix(env_name) / "bin" / binary


def create_env(manifest_path: Path, name: str) -> None:
    """Crea/reconstruye un entorno a partir de un manifiesto YAML."""
    ensure_pga_root()
    mm = find_micromamba()
    prefix = env_prefix(name)
    cmd = [mm, "create", "-y", "-p", str(prefix), "-f", str(manifest_path)]
    subprocess.run(cmd, check=True)


def install_package_editable(env_name: str, repo_root: Path) -> None:
    """Instala el paquete publicgenomicagent en modo editable dentro del entorno.

    Necesario porque micromamba create borra cualquier paquete pip previo
    que no esté declarado en el manifiesto.
    """
    pip = bin_path(env_name, "pip")
    if not pip.exists():
        raise FileNotFoundError(f"pip no encontrado en {env_name}: {pip}")
    subprocess.run([str(pip), "install", "-e", str(repo_root)], check=True)
