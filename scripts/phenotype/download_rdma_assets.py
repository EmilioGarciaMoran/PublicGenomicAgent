#!/usr/bin/env python3
"""Descarga los vector stores de RDMA desde el release de GitHub.

- Descarga EmbeddedDocs.zip si no está en cache.
- Verifica el SHA-256 del ZIP.
- Descomprime en ~/.pga/cache/rdma/.
- Verifica los SHA-256 de cada fichero individual.

Idempotente: si los ficheros ya están en cache con los hashes correctos,
no descarga nada.
"""
from __future__ import annotations

import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path


# Configuración
RELEASE_URL = (
    "https://github.com/EmilioGarciaMoran/PublicGenomicAgent/releases/download/"
    "v0.0.1-rdma-data/EmbeddedDocs.zip"
)
ZIP_SHA256 = "8980ac2c0b0b1f72d2b6e4926bda16998fd8198fe6bc182c541f1147ab17a642"

# Hashes de los ficheros individuales (verificados dentro del ZIP)
EXPECTED_FILES = {
    "vector_stores/G2GHPO_metadata_medembed.npy":
        "23acfecb52431348d73740bfd5df9a0e2d638c9cdb7f85f5a82d1fefe306debc",
    "vector_stores/rd_orpha_medembed.npy":
        "9e852750ed77e8e87a6c0e54906a91d45fd2ff91e4eee7d64863f30dd4ae9d85",
    "tools/abbreviations_medembed_sm.npy":
        "3fb49329528524d42693398a2b5e56b3d8e81b0116ff02ea740e7dcd9ffe04cc",
    "tools/lab_tables_medembed_sm.npy":
        "359112096a57590ad80e9992eb0ebd9e74a34b0c49e9a7657b9efc8705fb204a",
}

CACHE_ROOT = Path.home() / ".pga" / "cache" / "rdma"


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _all_files_present() -> bool:
    for rel, expected in EXPECTED_FILES.items():
        p = CACHE_ROOT / rel
        if not p.exists():
            return False
        if sha256_of_file(p) != expected:
            return False
    return True


def _download(url: str, dest: Path) -> None:
    print(f"  Descargando {url}")
    print(f"  → {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=300) as r:
        total = int(r.headers.get("content-length") or 0)
        downloaded = 0
        chunk_size = 1024 * 1024  # 1 MB
        with tmp.open("wb") as f:
            while True:
                chunk = r.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = 100 * downloaded / total
                    print(
                        f"\r  {pct:5.1f}%  ({downloaded//1024//1024} MB / "
                        f"{total//1024//1024} MB)",
                        end="", flush=True,
                    )
    print()
    tmp.rename(dest)


def _extract_and_verify(zip_path: Path) -> None:
    """Descomprime y verifica cada fichero."""
    print(f"  Descomprimiendo {zip_path.name}")
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path) as z:
        # El ZIP tiene una raíz EmbeddedDocs/
        for member in z.namelist():
            if not member.startswith("EmbeddedDocs/"):
                continue
            if member.endswith("/"):
                continue

            rel = member[len("EmbeddedDocs/"):]
            dest = CACHE_ROOT / rel
            dest.parent.mkdir(parents=True, exist_ok=True)

            with z.open(member) as src, dest.open("wb") as dst:
                dst.write(src.read())

    print("  Verificando hashes de ficheros")
    for rel, expected in EXPECTED_FILES.items():
        p = CACHE_ROOT / rel
        actual = sha256_of_file(p)
        if actual != expected:
            raise RuntimeError(
                f"Hash mismatch en {rel}\n"
                f"  esperado: {expected}\n"
                f"  obtenido: {actual}"
            )
        print(f"  [ok] {rel}")


def ensure_rdma_assets(verbose: bool = True) -> Path:
    """Garantiza que los vector stores están en cache. Devuelve CACHE_ROOT."""
    if _all_files_present():
        if verbose:
            print(f"  Vector stores ya presentes en {CACHE_ROOT}")
        return CACHE_ROOT

    zip_path = CACHE_ROOT / "_download" / "EmbeddedDocs.zip"

    # Si el ZIP ya existe con hash correcto, saltamos la descarga
    if zip_path.exists() and sha256_of_file(zip_path) == ZIP_SHA256:
        if verbose:
            print(f"  ZIP ya descargado: {zip_path}")
    else:
        if verbose:
            print(f"  Descargando vector stores (132 MB)...")
        _download(RELEASE_URL, zip_path)

        if verbose:
            print(f"  Verificando SHA-256 del ZIP")
        actual = sha256_of_file(zip_path)
        if actual != ZIP_SHA256:
            zip_path.unlink()
            raise RuntimeError(
                f"SHA-256 del ZIP no coincide.\n"
                f"  esperado: {ZIP_SHA256}\n"
                f"  obtenido: {actual}"
            )

    _extract_and_verify(zip_path)
    return CACHE_ROOT


def main() -> int:
    print("Descargando vector stores de RDMA")
    try:
        ensure_rdma_assets()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(f"OK: vector stores en {CACHE_ROOT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
