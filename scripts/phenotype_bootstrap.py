#!/usr/bin/env python3
"""Bootstrap de la capa de fenotipo de PublicGenomicAgent.

Descarga y configura:
1. Vector stores de RDMA (~300 MB)
2. LIRICAL (JAR + datos, ~450 MB)
3. Opcionalmente, el modelo Mistral 24B (~14 GB, requiere GPU)

Idempotente: ejecuciones sucesivas no hacen nada si todo está en cache.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Permitir importar los scripts hermanos
sys.path.insert(0, str(Path(__file__).parent))

from phenotype.download_rdma_assets import ensure_rdma_assets
from phenotype.install_lirical import ensure_lirical


def bootstrap(
    only_rdma: bool = False,
    only_lirical: bool = False,
    with_model: bool = False,
) -> int:
    print("=" * 60)
    print("  PublicGenomicAgent — phenotype bootstrap")
    print("=" * 60)

    rdma_ok = lirical_ok = model_ok = True

    if not only_lirical:
        print()
        print("[1/3] Vector stores de RDMA")
        try:
            ensure_rdma_assets()
            rdma_ok = True
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            rdma_ok = False

    if not only_rdma:
        print()
        print("[2/3] LIRICAL")
        try:
            ensure_lirical()
            lirical_ok = True
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            lirical_ok = False

    if with_model and not only_rdma and not only_lirical:
        print()
        print("[3/3] Modelo Mistral 24B (opcional)")
        print("  NOTA: la descarga del modelo no está implementada todavía.")
        print("  Ejecuta manualmente:")
        print("    huggingface-cli download mistralai/Mistral-Small-24B-Instruct-2501 \\")
        print("      --local-dir ~/.pga/cache/models/Mistral-Small-24B-Instruct-2501")

    print()
    print("=" * 60)
    if rdma_ok and lirical_ok:
        print("  Bootstrap completado correctamente")
        return 0
    else:
        print("  Bootstrap con errores")
        return 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Bootstrap de la capa de fenotipo de PublicGenomicAgent"
    )
    parser.add_argument(
        "--only-rdma", action="store_true",
        help="Solo descargar los vector stores de RDMA",
    )
    parser.add_argument(
        "--only-lirical", action="store_true",
        help="Solo instalar LIRICAL",
    )
    parser.add_argument(
        "--with-model", action="store_true",
        help="Descargar también Mistral 24B (~14 GB, requiere GPU)",
    )
    args = parser.parse_args()

    if args.only_rdma and args.only_lirical:
        print("ERROR: --only-rdma y --only-lirical son mutuamente excluyentes",
              file=sys.stderr)
        return 2

    return bootstrap(
        only_rdma=args.only_rdma,
        only_lirical=args.only_lirical,
        with_model=args.with_model,
    )


if __name__ == "__main__":
    sys.exit(main())
