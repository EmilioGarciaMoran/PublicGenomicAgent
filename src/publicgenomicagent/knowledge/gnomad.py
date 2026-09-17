from __future__ import annotations

import json
import urllib.request
from pathlib import Path

from ..env.paths import CACHE_DIR

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_CACHE = CACHE_DIR / "gnomad"

REGION_QUERY = """
query($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      pos
      ref
      alt
      exome { populations { id ac an } }
      genome { populations { id ac an } }
    }
  }
}
"""


def _cache_path(chrom: str, start: int, stop: int) -> Path:
    GNOMAD_CACHE.mkdir(parents=True, exist_ok=True)
    return GNOMAD_CACHE / f"{chrom}_{start}_{stop}.json"


def fetch_region(chrom: str, start: int, stop: int, refresh: bool = False) -> list[dict]:
    """Devuelve las variantes gnomAD v4 en la región. Cachea en disco."""
    cache = _cache_path(chrom, start, stop)
    if cache.exists() and not refresh:
        return json.loads(cache.read_text())

    body = json.dumps({
        "query": REGION_QUERY,
        "variables": {"chrom": chrom, "start": start, "stop": stop},
    }).encode()
    req = urllib.request.Request(
        GNOMAD_URL,
        data=body,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        payload = json.loads(r.read())

    variants = (payload.get("data", {})
                       .get("region", {})
                       .get("variants", [])) or []
    cache.write_text(json.dumps(variants))
    return variants


def mid_freq(variant: dict) -> tuple[int, int] | None:
    """Devuelve (AC, AN) para la población MID. Prefiere exome, luego genome.

    Retorna None si no hay datos MID o AN=0.
    """
    for source in ("exome", "genome"):
        block = variant.get(source)
        if not block:
            continue
        pops = block.get("populations") or []
        mid = next((p for p in pops if p["id"] == "mid"), None)
        if mid and mid["an"] > 0:
            return mid["ac"], mid["an"]
    return None


def filter_mid_variants(variants: list[dict], min_ac: int = 1) -> list[dict]:
    """Devuelve variantes con AC_MID >= min_ac."""
    out = []
    for v in variants:
        mid = mid_freq(v)
        if mid and mid[0] >= min_ac:
            out.append(v)
    return out
