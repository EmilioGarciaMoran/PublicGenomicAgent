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


# ---------------------------------------------------------------------------
# Frecuencias para el filtro de variantes raras (E1)
# ---------------------------------------------------------------------------
#
# Consulta con totales por bloque (exome/genome). Si el servidor rechazara
# esos campos, se cae a la consulta por poblaciones y el total se calcula
# sumando los grupos de ascendencia (sin los subconjuntos por sexo).

import time
import urllib.error

AF_QUERY = REGION_QUERY.replace(
    "exome { populations { id ac an } }", "exome { ac an populations { id ac an } }"
).replace(
    "genome { populations { id ac an } }", "genome { ac an populations { id ac an } }"
)
AF_CACHE = CACHE_DIR / "gnomad_af"
MAX_WINDOW = 50_000          # ventanas de consulta acotadas


def _post(query: str, chrom: str, start: int, stop: int, retries: int = 4) -> dict:
    body = json.dumps({"query": query,
                       "variables": {"chrom": chrom, "start": start, "stop": stop}}).encode()
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(GNOMAD_URL, data=body,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:     # 429: límite de peticiones
            last = e
            if e.code not in (429, 502, 503, 504):
                raise
        except urllib.error.URLError as e:
            last = e
        time.sleep(2 ** attempt * 5)
    raise RuntimeError(f"gnomAD no responde ({chrom}:{start}-{stop}): {last}")


def fetch_region_af(chrom: str, start: int, stop: int, refresh: bool = False) -> list[dict]:
    """Variantes gnomAD v4 con totales y poblaciones, en ventanas de <= 50 kb, con caché.

    `chrom` sin prefijo "chr". Tras la primera consulta funciona sin red.
    """
    chrom = chrom.removeprefix("chr")
    AF_CACHE.mkdir(parents=True, exist_ok=True)
    out: list[dict] = []
    for a in range(start, stop + 1, MAX_WINDOW):
        b = min(stop, a + MAX_WINDOW - 1)
        cache = AF_CACHE / f"{chrom}_{a}_{b}.json"
        if cache.exists() and not refresh:
            out += json.loads(cache.read_text())
            continue
        payload = _post(AF_QUERY, chrom, a, b)
        if payload.get("errors"):                 # campos no aceptados: consulta básica
            payload = _post(REGION_QUERY, chrom, a, b)
        if payload.get("errors"):
            raise RuntimeError(f"gnomAD devolvió errores: {payload['errors'][:2]}")
        variants = (payload.get("data") or {}).get("region", {}).get("variants") or []
        cache.write_text(json.dumps(variants))
        out += variants
    return out


SEX_SPLITS = {"XX", "XY"}


def _block_counts(block: dict | None) -> tuple[int, int, int, int, str] | None:
    """(AC, AN, AC_mid, AN_mid, método) de un bloque exome o genome."""
    if not block:
        return None
    pops = block.get("populations") or []
    mid = next((p for p in pops if p.get("id") == "mid"), None)
    ac_mid, an_mid = (mid["ac"], mid["an"]) if mid else (0, 0)
    if block.get("an") is not None and block.get("ac") is not None:
        return block["ac"], block["an"], ac_mid, an_mid, "total"
    groups = [p for p in pops if "_" not in p["id"] and p["id"] not in SEX_SPLITS]
    return (sum(p["ac"] for p in groups), sum(p["an"] for p in groups),
            ac_mid, an_mid, "suma_poblaciones")


def variant_frequencies(v: dict) -> dict:
    """AF global y AF de Oriente Medio (MID), sumando exoma y genoma de gnomAD v4."""
    ac = an = ac_m = an_m = 0
    methods = set()
    for key in ("exome", "genome"):
        c = _block_counts(v.get(key))
        if c:
            ac += c[0]; an += c[1]; ac_m += c[2]; an_m += c[3]; methods.add(c[4])
    return {
        "ac": ac, "an": an, "af": ac / an if an else None,
        "ac_mid": ac_m, "an_mid": an_m, "af_mid": ac_m / an_m if an_m else None,
        "method": "+".join(sorted(methods)) or "sin_datos",
    }
