from __future__ import annotations

import json
import urllib.request
from pathlib import Path

from ..env.paths import CACHE_DIR

UCSC_API = "https://api.genome.ucsc.edu/getData/sequence"
REF_CACHE = CACHE_DIR / "reference"


def _cache_path(genome: str, chrom: str, start: int, end: int) -> Path:
    REF_CACHE.mkdir(parents=True, exist_ok=True)
    return REF_CACHE / f"{genome}_{chrom}_{start}_{end}.fa"


def fetch_sequence(
    chrom: str,
    start: int,
    end: int,
    genome: str = "hg38",
    refresh: bool = False,
) -> Path:
    """Descarga la secuencia genómica de UCSC para una región y la cachea.

    UCSC usa 0-based half-open. La secuencia se uppercasea para
    evitar soft-masking (minúsculas), que rompe comparaciones de REF en VCFs.
    """
    cache = _cache_path(genome, chrom, start, end)
    if cache.exists() and not refresh:
        return cache

    api_start = start - 1
    url = f"{UCSC_API}?genome={genome};chrom={chrom};start={api_start};end={end}"
    with urllib.request.urlopen(url, timeout=60) as r:
        payload = json.loads(r.read())

    dna = payload.get("dna", "")
    if not dna:
        raise RuntimeError(f"UCSC no devolvió secuencia para {chrom}:{start}-{end}")

    dna = dna.upper()

    with cache.open("w") as f:
        f.write(f">{chrom}\n")
        for i in range(0, len(dna), 60):
            f.write(dna[i:i + 60] + "\n")

    with open(str(cache) + ".fai", "w") as f:
        header_len = len(chrom) + 2
        f.write(f"{chrom}\t{len(dna)}\t{header_len}\t60\t61\n")

    return cache
