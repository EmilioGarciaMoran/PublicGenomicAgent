"""Modelo de pedigree mínimo y construcción del fixture estándar.

Estructura del fixture (10 individuos, 3 generaciones):

       GF1--GM1     GF2--GM2       (fundadores)
          |            |
          F1 --------- M1          (padres)
             |  |  |  |
            C1 C2 C3 C4            (hijos: C2 y C3 afectados)

- C2: SV de novo (no presente en padres)
- C3: SV homocigoto recesivo (padres heterocigotos portadores)
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Individual:
    iid: str
    sex: str          # "1" = male, "2" = female, "0" = unknown
    affected: bool = False
    father: str | None = None
    mother: str | None = None


@dataclass
class Pedigree:
    individuals: dict[str, Individual] = field(default_factory=dict)

    def add(self, ind: Individual) -> None:
        self.individuals[ind.iid] = ind

    def children_of(self, parent_iid: str) -> list[str]:
        return [i.iid for i in self.individuals.values() if parent_iid in (i.father, i.mother)]

    def order_for_vcf(self) -> list[str]:
        """Fundadores primero, luego descendencia, orden natural VCF."""
        return ["GF1", "GM1", "GF2", "GM2", "F1", "M1", "C1", "C2", "C3", "C4"]


def build_standard_pedigree() -> Pedigree:
    p = Pedigree()
    # Fundadores
    p.add(Individual("GF1", "1"))
    p.add(Individual("GM1", "2"))
    p.add(Individual("GF2", "1"))
    p.add(Individual("GM2", "2"))
    # Padres
    p.add(Individual("F1", "1", father="GF1", mother="GM1"))
    p.add(Individual("M1", "2", father="GF2", mother="GM2"))
    # Hijos
    p.add(Individual("C1", "1", father="F1", mother="M1"))
    p.add(Individual("C2", "2", affected=True, father="F1", mother="M1"))
    p.add(Individual("C3", "1", affected=True, father="F1", mother="M1"))
    p.add(Individual("C4", "2", father="F1", mother="M1"))
    return p
