"""Calibración de la verificación familiar por submuestreo de sitios.

A partir de un VCF de genotipos en un panel de sitios (salida de
`verify_relatedness`), simula paneles de distintos tamaños tomando
subconjuntos aleatorios de los sitios informativos y mide, para cada
escenario declarado (un PED), con qué frecuencia la verificación da
cada estado.

Escenarios típicos con GIAB:
  - trío verdadero (HG003 padre, HG004 madre, HG002 hijo): debe pasar;
    la tasa de "failed" es la tasa de falsas alarmas.
  - padre falso (HG005 declarado padre de HG002): debe fallar;
    la tasa de "failed" es la sensibilidad.

Además resume, con todos los sitios, IBS0/n de cada pareja
padre/madre-hijo verdadera (suelo de error de genotipado) y el
relatedness de las parejas no emparentadas (fondo poblacional).
"""
from __future__ import annotations

import random
from collections import Counter
from dataclasses import dataclass

from .base import IndividualSpec
from .relatedness import evaluate, pair_stats

Genotypes = dict[tuple[str, int], dict[str, int | None]]


@dataclass
class Scenario:
    name: str
    pedigree: list[IndividualSpec]
    expected: str            # estado correcto: "passed" o "failed"


def informative_keys(gts: Genotypes, samples: list[str]) -> list[tuple[str, int]]:
    """Sitios genotipados (DP suficiente) en todas las muestras indicadas."""
    return [k for k, row in gts.items() if all(row.get(s) is not None for s in samples)]


def calibrate(
    gts: Genotypes,
    samples: list[str],
    scenarios: list[Scenario],
    ks: list[int],
    reps: int = 500,
    seed: int = 1,
    **thresholds,
) -> dict:
    """Tabla escenario x tamaño de panel -> frecuencia de cada estado."""
    rng = random.Random(seed)
    out: dict = {"n_sites_total": len(gts), "scenarios": {}}
    for sc in scenarios:
        members = [i.sample for i in sc.pedigree if i.sample in samples]
        keys = informative_keys(gts, members)
        rows = []
        for k in ks:
            if k > len(keys):
                rows.append({"k": k, "skipped": f"solo {len(keys)} sitios informativos"})
                continue
            c: Counter = Counter()
            for _ in range(reps):
                sub = {key: gts[key] for key in rng.sample(keys, k)}
                status, *_ = evaluate(sc.pedigree, sub, members, **thresholds)
                c[status] += 1
            rows.append({
                "k": k,
                "reps": reps,
                "counts": dict(c),
                "correct_rate": c[sc.expected] / reps,
                "failed_rate": c["failed"] / reps,
                "not_verifiable_rate": c["not_verifiable"] / reps,
            })
        full_status, full_pairs, issues, warnings = evaluate(
            sc.pedigree, {k: gts[k] for k in keys}, members, **thresholds)
        out["scenarios"][sc.name] = {
            "expected": sc.expected,
            "members": members,
            "n_informative": len(keys),
            "all_sites_status": full_status,
            "all_sites_pairs": full_pairs,
            "by_panel_size": rows,
        }
    return out


def pairwise_summary(gts: Genotypes, samples: list[str]) -> list[dict]:
    """IBS0/n y relatedness de todas las parejas con todos los sitios comunes."""
    res = []
    for i, a in enumerate(samples):
        for b in samples[i + 1:]:
            keys = informative_keys(gts, [a, b])
            st = pair_stats([gts[k][a] for k in keys], [gts[k][b] for k in keys])
            st["ibs0_rate"] = st["ibs0"] / st["n"] if st["n"] else None
            res.append({"a": a, "b": b, **st})
    return res


def markdown_report(cal: dict, pairs: list[dict], meta: dict) -> str:
    """Informe legible: parejas con todos los sitios + tabla por tamaño de panel."""
    lines = ["# Calibración de la verificación familiar (GIAB)", ""]
    for k, v in meta.items():
        lines.append(f"- **{k}:** {v}")
    lines += ["", "## Parejas con todos los sitios informativos", "",
              "| Pareja | n | IBS0 | IBS0/n | relatedness |", "|---|---|---|---|---|"]
    for p in pairs:
        r = "—" if p["relatedness"] is None else f"{p['relatedness']:.3f}"
        rate = "—" if p["ibs0_rate"] is None else f"{p['ibs0_rate']:.4f}"
        lines.append(f"| {p['a']}–{p['b']} | {p['n']} | {p['ibs0']} | {rate} | {r} |")
    for name, sc in cal["scenarios"].items():
        lines += ["", f"## Escenario: {name} (correcto = {sc['expected']})", "",
                  f"Sitios informativos: {sc['n_informative']}; con todos: "
                  f"**{sc['all_sites_status']}**.", "",
                  "| Sitios del panel | passed | warning | failed | not_verifiable | correcto |",
                  "|---|---|---|---|---|---|"]
        for row in sc["by_panel_size"]:
            if "skipped" in row:
                lines.append(f"| {row['k']} | {row['skipped']} | | | | |")
                continue
            c, n = row["counts"], row["reps"]
            rate = {s: c.get(s, 0) / n for s in ("passed", "warning", "failed", "not_verifiable")}
            lines.append(f"| {row['k']} | {rate['passed']:.3f} | {rate['warning']:.3f} | "
                         f"{rate['failed']:.3f} | {rate['not_verifiable']:.3f} | "
                         f"**{row['correct_rate']:.3f}** |")
    lines.append("")
    return "\n".join(lines)
