from __future__ import annotations

import json
from pathlib import Path

from .base import IndividualSpec


def load_pedigree_from_fam(fam_path: Path) -> list[IndividualSpec]:
    """Carga un .fam de PLINK y devuelve la lista de IndividualSpec.

    Formato .fam: FID IID PID MID sex phenotype
    - sex: 1=male, 2=female, 0=unknown
    - phenotype: 1=unaffected, 2=affected, 0/-9=missing
    """
    fam = Path(fam_path).expanduser().resolve()
    if not fam.exists():
        raise FileNotFoundError(f".fam no encontrado: {fam}")

    rows: list[tuple[str, str, str, str, str, str]] = []
    with fam.open() as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 6:
                continue
            rows.append(tuple(fields[:6]))

    specs: list[IndividualSpec] = []
    for fid, iid, pid, mid, sex, pheno in rows:
        sex_str = {"1": "M", "2": "F"}.get(sex, "U")
        affected = (pheno == "2")
        father = None if pid in ("0", "-9") else pid
        mother = None if mid in ("0", "-9") else mid
        specs.append(IndividualSpec(
            sample=iid,
            sex=sex_str,
            affected=affected,
            father=father,
            mother=mother,
        ))
    return specs


def load_pedigree_from_json(json_path: Path) -> list[IndividualSpec]:
    """Carga un pedigrí desde JSON estructurado.

    Formato esperado:
        {
          "individuals": [
            {"sample": "C2", "sex": "M", "affected": true,
             "father": "F1", "mother": "M1"},
            ...
          ]
        }
    """
    p = Path(json_path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"JSON de pedigrí no encontrado: {p}")

    data = json.loads(p.read_text())
    individuals = data.get("individuals", [])

    specs: list[IndividualSpec] = []
    for ind in individuals:
        specs.append(IndividualSpec(
            sample=ind["sample"],
            sex=ind.get("sex", "U"),
            affected=bool(ind.get("affected", False)),
            father=ind.get("father"),
            mother=ind.get("mother"),
        ))
    return specs


def find_proband(pedigree: list[IndividualSpec]) -> str:
    """Heurística: devuelve el sample del único afectado, si hay uno solo."""
    affected = [i for i in pedigree if i.affected]
    if len(affected) == 1:
        return affected[0].sample
    raise ValueError(
        f"No se puede inferir el probando: {len(affected)} afectados. "
        "Especifica --proband explícitamente."
    )


def get_parents(
    pedigree: list[IndividualSpec], sample: str,
) -> tuple[str | None, str | None]:
    """Devuelve (father_sample, mother_sample) para un sample dado."""
    for ind in pedigree:
        if ind.sample == sample:
            return ind.father, ind.mother
    raise ValueError(f"Sample '{sample}' no encontrado en el pedigrí")
