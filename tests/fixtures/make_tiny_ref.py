#!/usr/bin/env python3
"""Genera una referencia FASTA sintética de chr1:1-10000 con fondo 'A'.

Los reads del fixture make_tiny_bam.sh tienen bases distintas al fondo,
así que el variant calling encontrará variantes en las regiones:
- chr1:100-209   (reads r1-r5)
- chr1:500-609   (reads r6-r8)
- chr1:9000-9109 (reads r9-r10)

Uso: make_tiny_ref.py <out.fa>
"""
from __future__ import annotations

import sys
from pathlib import Path

LENGTH = 10000
LINE_WIDTH = 60


def main() -> None:
    if len(sys.argv) != 2:
        print("Uso: make_tiny_ref.py <out.fa>", file=sys.stderr)
        sys.exit(2)

    out = Path(sys.argv[1])
    seq = "A" * LENGTH

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        f.write(">chr1\n")
        for i in range(0, len(seq), LINE_WIDTH):
            f.write(seq[i:i + LINE_WIDTH] + "\n")

    # .fai: nombre, longitud, offset, linebases, linewidth
    fai = Path(str(out) + ".fai")
    with fai.open("w") as f:
        # header ">chr1\n" = 6 bytes, offset de la primera base = 6
        f.write(f"chr1\t{LENGTH}\t6\t{LINE_WIDTH}\t{LINE_WIDTH + 1}\n")

    print(f"OK: {out}")


if __name__ == "__main__":
    main()
