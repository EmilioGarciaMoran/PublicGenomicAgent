"""Catálogo de variantes estructurales para el fixture.

Las variantes son sintéticas pero plausibles por tipo y tamaño:
- DEL 150 pb (microdeleción exónica)
- DEL 2 kb (exón completo)
- DEL 8 kb (multiexónica)
- INS 300 pb (tipo Alu)
- INS 6 kb (tipo LINE-1)
- DUP 5 kb (tandem)
- INV 3 kb (inversión)
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SV:
    id: str
    chrom: str
    pos: int           # posición 1-based del breakpoint
    svtype: str        # DEL | INS | DUP | INV
    end: int           # para DEL/DUP/INV: fin del evento; para INS: pos (mismo)
    length: int        # tamaño
    in_children: tuple[str, ...] = ()
    inheritance: str = "common"  # common | rare | de_novo


def default_catalog(roi_chrom: str, roi_start: int) -> list[SV]:
    """7 SVs distribuidos a lo largo del ROI NPHP1 (~25 kb)."""
    base = roi_start
    return [
        SV("SV001", roi_chrom, base + 2000, "DEL", base + 2150, 150),
        SV("SV002", roi_chrom, base + 5000, "DEL", base + 7000, 2000),
        SV("SV003", roi_chrom, base + 8000, "DEL", base + 16000, 8000),
        SV("SV004", roi_chrom, base + 18000, "INS", base + 18000, 300),
        SV("SV005", roi_chrom, base + 20000, "INS", base + 20000, 6000),
        SV("SV006", roi_chrom, base + 12000, "DUP", base + 17000, 5000),
        SV("SV007", roi_chrom, base + 3000, "INV", base + 6000, 3000),
    ]
