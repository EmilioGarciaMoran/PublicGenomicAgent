from __future__ import annotations

import gzip
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from ..env.micromamba import bin_path
from ..env.runtime import ToolRuntime
from .base import (
    MendelianFilterInput,
    MendelianFilterOutput,
    PlinkValidateInput,
    PlinkValidateOutput,
)


# ---------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------

def _bcftools() -> str:
    p = bin_path("pga-hts", "bcftools")
    if not p.exists():
        raise FileNotFoundError(f"bcftools no encontrado: {p}")
    return str(p)


def _plink() -> str:
    p = bin_path("pga-mendelian", "plink")
    if not p.exists():
        raise FileNotFoundError(f"plink no encontrado: {p}")
    return str(p)


@dataclass
class SampleData:
    name: str
    gt: str          # "0/1", "1/1", "0/0", "./.", etc.
    dp: int | None


def _parse_gt(gt_field: str) -> tuple[str, int | None]:
    """Extrae GT y DP de un campo FORMAT/muestra."""
    if ":" not in gt_field:
        return gt_field, None
    parts = gt_field.split(":")
    gt = parts[0]
    # No tenemos acceso a los índices del FORMAT aquí; se resuelve en el caller
    return gt, None


def _sample_gt(fields: list[str], samples: list[str], sample: str,
               fmt_keys: list[str]) -> tuple[str, int | None]:
    """Devuelve (GT, DP) para un sample concreto."""
    idx = samples.index(sample)
    sample_field = fields[9 + idx].split(":")
    gt = sample_field[0] if sample_field else "./."
    dp = None
    if "DP" in fmt_keys:
        dp_idx = fmt_keys.index("DP")
        if dp_idx < len(sample_field) and sample_field[dp_idx].isdigit():
            dp = int(sample_field[dp_idx])
    return gt, dp


def _iter_vcf(vcf: Path):
    """Iterador sobre líneas de un VCF, gz o plano."""
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt") as f:
        for line in f:
            yield line


def _write_vcf_header(out_f, header_lines: list[str]) -> None:
    for line in header_lines:
        if line.startswith("##"):
            out_f.write(line)
    for line in header_lines:
        if line.startswith("#CHROM"):
            out_f.write(line)
            break
