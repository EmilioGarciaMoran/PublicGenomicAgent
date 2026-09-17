"""Simulación de genotipos: fundadores por Hardy-Weinberg,
descendientes por segregación mendeliana.
"""
from __future__ import annotations

import random

from pedigree_model import Pedigree


def founder_genotype(af: float, rng: random.Random) -> tuple[int, int]:
    """Muestrea genotipo 0/0, 0/1, 1/1 según AF (asume HWE)."""
    a1 = 1 if rng.random() < af else 0
    a2 = 1 if rng.random() < af else 0
    return (a1, a2)


def transmit(parent_gt: tuple[int, int], rng: random.Random) -> int:
    """Transmite un alelo al azar."""
    return parent_gt[rng.randint(0, 1)]


def simulate_snv_genotypes(
    pedigree: Pedigree,
    variants: list[dict],
    rng: random.Random,
) -> dict[str, dict[str, tuple[int, int]]]:
    """Devuelve {iid: {variant_id: (a1, a2)}}."""
    result: dict[str, dict[str, tuple[int, int]]] = {
        iid: {} for iid in pedigree.individuals
    }

    for v in variants:
        vid = v["variant_id"]
        af = v.get("_af_mid", 0.0)
        # Fundadores
        for iid, ind in pedigree.individuals.items():
            if ind.father is None and ind.mother is None:
                result[iid][vid] = founder_genotype(af, rng)
        # Descendientes en orden topológico
        for iid in pedigree.order_for_vcf():
            ind = pedigree.individuals[iid]
            if ind.father is None:
                continue
            f_gt = result[ind.father].get(vid, (0, 0))
            m_gt = result[ind.mother].get(vid, (0, 0))
            result[iid][vid] = (transmit(f_gt, rng), transmit(m_gt, rng))
    return result


def simulate_sv_genotypes(
    pedigree: Pedigree,
    catalog,
    rng: random.Random,
) -> dict[str, dict[str, tuple[int, int]]]:
    """Genotipos para SVs.

    - SV002: de novo en C2 (padres 0/0).
    - SV005: recesivo, C3 homocigoto, F1 y M1 heterocigotos.
    - Resto: al azar en un fundador.
    """
    result: dict[str, dict[str, tuple[int, int]]] = {
        iid: {} for iid in pedigree.individuals
    }
    for sv in catalog:
        for iid in pedigree.individuals:
            result[iid][sv.id] = (0, 0)

        if sv.id == "SV002":
            result["C2"][sv.id] = (0, 1)
        elif sv.id == "SV005":
            result["F1"][sv.id] = (0, 1)
            result["M1"][sv.id] = (0, 1)
            result["C3"][sv.id] = (1, 1)
        else:
            if rng.random() < 0.4:
                founder = rng.choice(["GF1", "GM1", "GF2", "GM2"])
                result[founder][sv.id] = (0, 1)
    return result
