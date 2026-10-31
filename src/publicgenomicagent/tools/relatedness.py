"""verify_relatedness: verificación familiar con IBS calculado en local.

Sin dependencias externas nuevas: genotipa un conjunto fijo de sitios
(SNV comunes) con `bcftools` (ya en pga-hts) en los BAMs del caso y
compara, por cada pareja de individuos secuenciados, lo observado con
lo declarado en el PED.

Métricas por pareja (definiciones de somalier, Pedersen et al. 2020,
reimplementadas aquí):
  n            sitios genotipados en ambos (DP >= min_dp)
  ibs0         uno 0/0 y el otro 1/1 (homocigotos opuestos)
  ibs2         mismo genotipo
  shared_hets  ambos 0/1
  relatedness  (shared_hets - 2*ibs0) / min(hets_i, hets_j)
               ~0 no emparentados, ~0.5 padre-hijo y hermanos, ~1 duplicado

Comprobaciones (calibradas con exomas GIAB, ver
scripts/calibrate_relatedness_giab.py):
  - padre/madre-hijo: prueba binomial sobre IBS0. Con un progenitor
    verdadero, IBS0 solo aparece por error de genotipado (GIAB: 0.03% y
    0.2% en las dos parejas reales); se toma genotype_error = 0.5%. El
    umbral es el menor recuento con P(IBS0 >= umbral | progenitor) < alpha
    (1e-3). IBS0 >= umbral -> FALLO. Es la prueba más robusta: un
    progenitor y su hijo no pueden ser homocigotos opuestos, haya
    consanguinidad o no.
  - potencia: si con n sitios la probabilidad de detectar a un
    no-progenitor DE LA MISMA POBLACIÓN (tasa de IBS0 de referencia 8.7%,
    la del matrimonio askenazí HG003-HG004, población endogámica) es
    < min_power (99.9%), la pareja es NO VERIFICABLE. En la práctica hacen
    falta ~200 sitios informativos. La regla fija anterior (IBS0/n > 1%)
    daba 1-4% de falsas alarmas con 100-150 sitios en el trío GIAB real.
  - relatedness >= dup_relatedness entre individuos distintos -> FALLO
    (duplicado o intercambio; el PED no expresa gemelos MZ).
  - parentesco esperado >= 0.25 (no padre/madre-hijo) y observado
    < esperado - deficit_tol -> FALLO (p.ej. hermanos que no lo son).
  - esperado < 0.1 y observado > undeclared_relatedness -> AVISO
    (parentesco no declarado; fondo askenazí observado: 0.08).
  - relatedness solo se estima con >= min_hets heterocigotos (ruido).
  - n < min_sites -> NO VERIFICABLE (suelo absoluto; no equivale a aprobada).
    Los sitios deben ser comunes (MAF alta): con MAF baja IBS0 es raro.

Una incoherencia se describe como incompatibilidad con el pedigrí
declarado, sin inferir causa (ver docs/mendelian_segregation.md, §7).
"""
from __future__ import annotations

import gzip
import subprocess
from functools import lru_cache
from itertools import combinations
from math import exp, fsum, lgamma, log, log1p
from pathlib import Path

from ..env.micromamba import bin_path
from ..env.runtime import ToolRuntime
from .base import (
    IndividualSpec,
    VerifyRelatednessInput,
    VerifyRelatednessOutput,
)
from .pedigree import kinship_function, relationship

Site = tuple[str, int, str, str]          # chrom, pos (1-based), ref, alt


# ---------------------------------------------------------------------------
# Sitios
# ---------------------------------------------------------------------------

def _open_text(p: Path):
    p = Path(p)
    return gzip.open(p, "rt") if p.suffix == ".gz" else p.open()


def load_sites(path: Path) -> list[Site]:
    """VCF(.gz) o TSV de 4 columnas: CHROM POS REF ALT. Solo SNV bialélicos."""
    sites: list[Site] = []
    with _open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) >= 5:                      # VCF: CHROM POS ID REF ALT ...
                chrom, pos, ref, alt = c[0], c[1], c[3], c[4]
            elif len(c) == 4:                    # TSV: CHROM POS REF ALT
                chrom, pos, ref, alt = c
            else:
                continue
            if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
                sites.append((chrom, int(pos), ref, alt))
    return sites


# ---------------------------------------------------------------------------
# Genotipos
# ---------------------------------------------------------------------------

def parse_genotypes(
    vcf_lines, sites: list[Site], min_dp: int,
) -> tuple[list[str], dict[tuple[str, int], dict[str, int | None]]]:
    """Dosis del alelo ALT del sitio (0/1/2) por muestra; None si no fiable.

    Descarta llamadas con DP < min_dp, alelos distintos del sitio o GT ausente.
    """
    want = {(c, p): (r, a) for c, p, r, a in sites}
    samples: list[str] = []
    out: dict[tuple[str, int], dict[str, int | None]] = {}
    for line in vcf_lines:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = f[9:]
            continue
        key = (f[0], int(f[1]))
        if key not in want:
            continue
        sref, salt = want[key]
        if f[3] != sref:
            continue
        alleles = [f[3]] + ([] if f[4] in (".", "<*>") else
                            [a for a in f[4].split(",") if a != "<*>"])
        fmt = f[8].split(":")
        gi = fmt.index("GT") if "GT" in fmt else None
        di = fmt.index("DP") if "DP" in fmt else None
        row: dict[str, int | None] = {}
        for s, val in zip(samples, f[9:]):
            v = val.split(":")
            dose: int | None = None
            gt = v[gi] if gi is not None and gi < len(v) else "."
            dp_ok = True
            if di is not None and di < len(v) and v[di] not in (".", ""):
                dp_ok = int(v[di]) >= min_dp
            elif di is not None:
                dp_ok = False
            idx = gt.replace("|", "/").split("/")
            if dp_ok and len(idx) == 2 and all(x.isdigit() for x in idx):
                al = [alleles[int(x)] if int(x) < len(alleles) else None for x in idx]
                if all(a in (sref, salt) for a in al):
                    dose = sum(1 for a in al if a == salt)
            row[s] = dose
        out[key] = row
    return samples, out


# ---------------------------------------------------------------------------
# Estadísticos por pareja
# ---------------------------------------------------------------------------

def pair_stats(gi: list[int | None], gj: list[int | None], min_hets: int = 10) -> dict:
    n = ibs0 = ibs2 = shared = hi = hj = 0
    for a, b in zip(gi, gj):
        if a is None or b is None:
            continue
        n += 1
        if a == b:
            ibs2 += 1
        if {a, b} == {0, 2}:
            ibs0 += 1
        if a == 1:
            hi += 1
        if b == 1:
            hj += 1
        if a == 1 and b == 1:
            shared += 1
    m = min(hi, hj)
    rel = (shared - 2 * ibs0) / m if m >= min_hets else None
    return {"n": n, "ibs0": ibs0, "ibs2": ibs2, "shared_hets": shared,
            "hets_a": hi, "hets_b": hj, "relatedness": rel}


def _log_pmf(k: int, n: int, p: float) -> float:
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
            + k * log(p) + (n - k) * log1p(-p))


@lru_cache(maxsize=65536)
def binom_sf(x: int, n: int, p: float) -> float:
    """P(X >= x) con X ~ Binomial(n, p), estable para n y x grandes (escala log)."""
    if x <= 0:
        return 1.0
    if x > n:
        return 0.0
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    if x > n * p:   # cola superior: sumar directamente k = x..n
        return min(1.0, fsum(exp(_log_pmf(k, n, p)) for k in range(x, n + 1)))
    return max(0.0, 1.0 - fsum(exp(_log_pmf(k, n, p)) for k in range(x)))


@lru_cache(maxsize=4096)
def ibs0_threshold(n: int, genotype_error: float, alpha: float) -> int:
    """Menor recuento de IBS0 improbable (< alpha) para un progenitor verdadero."""
    x = 0
    while binom_sf(x, n, genotype_error) >= alpha:
        x += 1
    return x


def evaluate(
    pedigree: list[IndividualSpec],
    genotypes: dict[tuple[str, int], dict[str, int | None]],
    samples: list[str],
    *,
    min_sites: int = 100,
    min_hets: int = 50,
    genotype_error: float = 0.005,
    alpha: float = 1e-3,
    nonparent_ibs0_rate: float = 0.087,
    min_power: float = 0.999,
    dup_relatedness: float = 0.85,
    deficit_tol: float = 0.2,
    undeclared_relatedness: float = 0.3,
) -> tuple[str, list[dict], list[str], list[str]]:
    """Compara observado con declarado. Devuelve (status, pares, issues, warnings)."""
    ped_ids = {i.sample for i in pedigree}
    seq = [s for s in samples if s in ped_ids]
    phi = kinship_function(pedigree)
    keys = list(genotypes)
    issues: list[str] = []
    warnings: list[str] = []
    pairs: list[dict] = []
    po_unverifiable = False
    for a, b in combinations(seq, 2):
        st = pair_stats([genotypes[k].get(a) for k in keys],
                        [genotypes[k].get(b) for k in keys], min_hets)
        rel_label = relationship(pedigree, a, b)
        exp = 2 * phi(a, b)
        r = st["relatedness"]
        status = "passed"
        notes: list[str] = []
        po = rel_label == "padre/madre-hijo"
        test: dict = {}
        if st["n"] < min_sites:
            status = "not_verifiable"
            notes.append(f"solo {st['n']} sitios informativos (< {min_sites})")
            if po:
                po_unverifiable = True
        else:
            if po:
                thr = ibs0_threshold(st["n"], genotype_error, alpha)
                power = binom_sf(thr, st["n"], nonparent_ibs0_rate)
                p_val = binom_sf(st["ibs0"], st["n"], genotype_error)
                test = {"ibs0_threshold": thr, "p_value": p_val, "power": round(power, 6)}
                if st["ibs0"] >= thr:
                    status = "failed"
                    notes.append(
                        f"incompatible con padre/madre-hijo declarado: IBS0 "
                        f"{st['ibs0']}/{st['n']} homocigotos opuestos (umbral {thr}, "
                        f"p={p_val:.1e} si fuera progenitor)")
                elif power < min_power:
                    status = "not_verifiable"
                    po_unverifiable = True
                    notes.append(
                        f"potencia insuficiente: con {st['n']} sitios se detectaría a un "
                        f"no-progenitor de la misma población con probabilidad "
                        f"{power:.3f} (< {min_power})")
            if r is not None and r >= dup_relatedness:
                status = "failed"
                notes.append(f"relatedness {r:.2f}: posible duplicado o intercambio de muestras")
            if r is not None and not po and exp >= 0.25 and r < exp - deficit_tol:
                status = "failed"
                notes.append(
                    f"parentesco observado {r:.2f} menor que el declarado "
                    f"({rel_label}, esperado {exp:.2f})")
            if (r is not None and exp < 0.1 and r > undeclared_relatedness
                    and status != "failed"):
                status = "warning"
                notes.append(f"parentesco no declarado: observado {r:.2f}, esperado {exp:.2f}")
            if r is None and status == "passed" and not po and exp >= 0.1:
                status = "not_verifiable"
                notes.append(f"pocos heterocigotos para estimar relatedness (< {min_hets})")
        pairs.append({"a": a, "b": b, "declared": rel_label,
                      "expected_relatedness": round(exp, 4), **st, **test,
                      "status": status, "notes": notes})
        for n_ in notes:
            (issues if status == "failed" else warnings).append(f"{a}–{b}: {n_}")

    if any(p["status"] == "failed" for p in pairs):
        overall = "failed"
    elif not pairs or all(p["status"] == "not_verifiable" for p in pairs) or po_unverifiable:
        overall = "not_verifiable"
        if not pairs:
            warnings.append("menos de dos individuos secuenciados del PED: nada que comparar")
    elif any(p["status"] == "warning" for p in pairs) or any(
            p["status"] == "not_verifiable" and p["expected_relatedness"] >= 0.1 for p in pairs):
        # Una pareja no emparentada que no se puede comprobar no aporta aviso.
        overall = "warning"
    else:
        overall = "passed"
    return overall, pairs, issues, warnings


# ---------------------------------------------------------------------------
# Tool
# ---------------------------------------------------------------------------

def genotype_command(bcftools: str, reference: Path, targets: Path,
                     bams: list[Path], out_vcf: Path) -> tuple[list[str], list[str]]:
    """mpileup en sitios fijos (sin -v: también 0/0) | call -m."""
    mp = [bcftools, "mpileup", "-B", "-q", "20", "-Q", "20", "-f", str(reference),
          "-T", str(targets), "-a", "FORMAT/DP,FORMAT/AD", *map(str, bams)]
    call = [bcftools, "call", "-m", "-Oz", "-o", str(out_vcf)]
    return mp, call


def verify_relatedness(runtime: ToolRuntime, inp: VerifyRelatednessInput) -> VerifyRelatednessOutput:
    out_dir = Path(inp.output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    sites = load_sites(inp.sites)
    targets = out_dir / "relatedness_targets.tsv"
    targets.write_text("".join(f"{c}\t{p}\n" for c, p, _, _ in sites))
    out_vcf = out_dir / "relatedness_sites.vcf.gz"

    bcftools = bin_path("pga-hts", "bcftools")
    if not bcftools.exists():
        raise FileNotFoundError(f"bcftools no encontrado: {bcftools}")
    mp, call = genotype_command(str(bcftools), Path(inp.reference_fasta), targets,
                                [Path(b) for _, b in inp.bams], out_vcf)
    p1 = subprocess.Popen(mp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(call, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p1.stdout:
        p1.stdout.close()
    _, e2 = p2.communicate()
    _, e1 = p1.communicate()
    if p2.returncode != 0 or p1.returncode not in (0, None):
        raise RuntimeError(f"bcftools falló: {(e1 or b'').decode()[-300:]} {(e2 or b'').decode()[-300:]}")

    with gzip.open(out_vcf, "rt") as f:
        samples, gts = parse_genotypes(f, sites, inp.min_dp)
    status, pairs, issues, warnings = evaluate(
        inp.pedigree, gts, samples,
        min_sites=inp.min_sites, min_hets=inp.min_hets,
        genotype_error=inp.genotype_error, alpha=inp.alpha,
        nonparent_ibs0_rate=inp.nonparent_ibs0_rate, min_power=inp.min_power,
        dup_relatedness=inp.dup_relatedness,
        deficit_tol=inp.deficit_tol, undeclared_relatedness=inp.undeclared_relatedness,
    )
    return VerifyRelatednessOutput(
        ok=True, tool="verify_relatedness",
        message=f"verificación familiar: {status} ({len(sites)} sitios, {len(pairs)} parejas)",
        status=status, n_sites_panel=len(sites), samples=samples, pairs=pairs,
        issues=issues, warnings=warnings, genotypes_vcf=out_vcf,
        outputs={"genotypes_vcf": str(out_vcf)},
    )
