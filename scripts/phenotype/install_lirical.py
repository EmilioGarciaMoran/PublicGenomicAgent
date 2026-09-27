#!/usr/bin/env python3
"""Descarga e instala LIRICAL en ~/.pga/cache/lirical/.

Pasos:
1. Descarga el ZIP de LIRICAL desde GitHub Releases.
2. Descomprime en ~/.pga/cache/lirical/.
3. Ejecuta `java -jar lirical-cli-*.jar download -d data` para bajar los datos.
4. Aplica el fix de duplicados en mim2gene_medgen y hgnc_complete_set.txt.

Idempotente: si todo está en su sitio, no hace nada.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path


LIRICAL_VERSION = "2.4.1"
RELEASE_URL = (
    f"https://github.com/TheJacksonLaboratory/LIRICAL/releases/download/"
    f"v{LIRICAL_VERSION}/lirical-cli-{LIRICAL_VERSION}-distribution.zip"
)

CACHE_ROOT = Path.home() / ".pga" / "cache" / "lirical"


def _find_jar() -> Path | None:
    if not CACHE_ROOT.exists():
        return None
    jars = list(CACHE_ROOT.rglob(f"lirical-cli-{LIRICAL_VERSION}.jar"))
    return jars[0] if jars else None


def _find_data_dir(jar: Path) -> Path | None:
    data = jar.parent / "data"
    return data if data.exists() else None


def _has_lirical_data(data_dir: Path) -> bool:
    """Verifica que los ficheros clave de datos existen."""
    required = [
        "mim2gene_medgen",
        "hgnc_complete_set.txt",
        "hp.json",
        "phenotype.hpoa",
    ]
    return all((data_dir / r).exists() for r in required)


def _check_java() -> str:
    java = shutil.which("java")
    if java is None:
        raise RuntimeError(
            "Java no encontrado. LIRICAL requiere Java 17+.\n"
            "Instálalo con: sudo apt install openjdk-17-jre"
        )
    return java


def _download(url: str, dest: Path) -> None:
    print(f"  Descargando {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=300) as r:
        total = int(r.headers.get("content-length") or 0)
        downloaded = 0
        with tmp.open("wb") as f:
            while True:
                chunk = r.read(1024 * 1024)
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


def _install_lirical_zip() -> Path:
    """Descarga y descomprime LIRICAL. Devuelve el path del JAR."""
    jar = _find_jar()
    if jar is not None:
        print(f"  LIRICAL ya instalado: {jar}")
        return jar

    zip_path = CACHE_ROOT / "_download" / f"lirical-cli-{LIRICAL_VERSION}.zip"
    if not zip_path.exists():
        _download(RELEASE_URL, zip_path)

    print(f"  Descomprimiendo {zip_path.name}")
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as z:
        z.extractall(CACHE_ROOT)

    jar = _find_jar()
    if jar is None:
        raise RuntimeError(
            f"No se encontró lirical-cli-{LIRICAL_VERSION}.jar tras descomprimir"
        )
    return jar


def _download_lirical_data(jar: Path) -> Path:
    """Ejecuta java -jar lirical download -d data."""
    data_dir = jar.parent / "data"
    if _has_lirical_data(data_dir):
        print(f"  Datos de LIRICAL ya presentes en {data_dir}")
        return data_dir

    java = _check_java()
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"  Descargando datos de LIRICAL (HPO, OMIM, Orphanet, ~400 MB)")
    result = subprocess.run(
        [java, "-jar", str(jar), "download", "-d", str(data_dir)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"LIRICAL download falló (exit {result.returncode}).\n"
            f"stderr:\n{result.stderr[-2000:]}"
        )
    return data_dir


def _apply_duplicate_fix(data_dir: Path) -> None:
    """Aplica el fix de duplicados en mim2gene_medgen y hgnc_complete_set.txt.

    LIRICAL 2.4.1 falla al arrancar con los datos actuales de HGNC/NCBI
    porque Phenol construye un mapa con Collectors.toMap y hay claves
    duplicadas. Filtramos manteniendo solo la primera aparición.
    """
    # Fix 1: mim2gene_medgen (columna 2 = GeneID)
    mim = data_dir / "mim2gene_medgen"
    if mim.exists():
        print(f"  Aplicando fix de duplicados en mim2gene_medgen")
        lines = mim.read_text().splitlines()
        seen: set[str] = set()
        out: list[str] = []
        for line in lines:
            if line.startswith("#"):
                out.append(line)
                continue
            fields = line.split("\t")
            if len(fields) < 2 or fields[1] == "-":
                out.append(line)
                continue
            if fields[1] in seen:
                continue
            seen.add(fields[1])
            out.append(line)
        mim.write_text("\n".join(out) + "\n")

    # Fix 2: hgnc_complete_set.txt (columna 19 = ensembl_gene_id)
    hgnc = data_dir / "hgnc_complete_set.txt"
    if hgnc.exists():
        print(f"  Aplicando fix de duplicados en hgnc_complete_set.txt")
        lines = hgnc.read_text().splitlines()
        seen_h: set[str] = set()
        out_h: list[str] = []
        for i, line in enumerate(lines):
            if i == 0:
                out_h.append(line)
                continue
            fields = line.split("\t")
            if len(fields) < 19 or not fields[18]:
                out_h.append(line)
                continue
            if fields[18] in seen_h:
                continue
            seen_h.add(fields[18])
            out_h.append(line)
        hgnc.write_text("\n".join(out_h) + "\n")


def ensure_lirical(verbose: bool = True) -> Path:
    """Garantiza que LIRICAL está instalado y con datos. Devuelve el JAR."""
    jar = _install_lirical_zip()
    data_dir = _download_lirical_data(jar)
    _apply_duplicate_fix(data_dir)
    return jar


def main() -> int:
    print("Instalando LIRICAL")
    try:
        jar = ensure_lirical()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(f"OK: LIRICAL en {jar}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
