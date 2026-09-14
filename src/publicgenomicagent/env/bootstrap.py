from __future__ import annotations

from pathlib import Path

from .manifests import read_lock, sha256_of_file, write_lock
from .micromamba import create_env, env_exists, install_package_editable
from .registry import Registry


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def bootstrap_env(registry: Registry, name: str, force: bool = False) -> str:
    if name not in registry.envs:
        raise KeyError(f"Entorno '{name}' no declarado en registry.yaml")

    spec = registry.envs[name]
    manifest = spec.manifest

    if not manifest.exists():
        raise FileNotFoundError(f"Manifiesto no encontrado: {manifest}")

    current_hash = sha256_of_file(manifest)
    lock = read_lock(name)
    exists = env_exists(name)

    if exists and lock and lock.get("manifest_hash") == current_hash and not force:
        return f"[=] {name} ya está actualizado ({current_hash[:12]})"

    # Caso: entorno existe pero no hay lock → solo registrar e instalar pkg si aplica.
    if exists and lock is None and not force:
        if spec.installs_self:
            install_package_editable(name, repo_root())
        write_lock(name, {
            "env": name,
            "manifest": str(manifest),
            "manifest_hash": current_hash,
        })
        return f"[+] {name} registrado (sin lock previo) — hash {current_hash[:12]}"

    # Caso: crear/reconstruir entorno desde manifiesto, luego instalar pkg si aplica.
    create_env(manifest, name)
    if spec.installs_self:
        install_package_editable(name, repo_root())

    write_lock(name, {
        "env": name,
        "manifest": str(manifest),
        "manifest_hash": current_hash,
    })

    action = "reconstruyendo (forzado)" if force else (
        "manifiesto cambió, reconstruyendo" if exists else "creando"
    )
    return f"[+] {name} {action} — hash {current_hash[:12]}"
