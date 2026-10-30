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


# ---------------------------------------------------------------------------
# PED: lectura, escritura y validación estructural
# ---------------------------------------------------------------------------
#
# El PED es la única fuente del pedigrí. Puede generarse desde texto o
# desde una foto, pero siempre como borrador que una persona confirma, y
# siempre pasa por `validate_pedigree` antes de usarse. La consanguinidad
# se expresa dibujando el bucle: los antepasados comunes aparecen como
# filas del PED (aunque no estén secuenciados).

from dataclasses import dataclass, field
from functools import lru_cache


@dataclass
class PedValidation:
    issues: list[str] = field(default_factory=list)     # bloquean
    warnings: list[str] = field(default_factory=list)   # se registran

    @property
    def ok(self) -> bool:
        return not self.issues


def parse_ped_text(text: str) -> tuple[list[str], list[IndividualSpec]]:
    """Parsea un PED (6 columnas: FID IID PID MID SEX PHENO).

    SEX: 1=M, 2=F, otro=U. PHENO: 2=afectado, 1=no afectado,
    0/-9/otro=desconocido (phenotype_known=False, affected=False).
    Ignora líneas vacías y comentarios '#'. Devuelve (family_ids, individuos).
    """
    fids: list[str] = []
    specs: list[IndividualSpec] = []
    for n, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 6:
            raise ValueError(f"PED línea {n}: se esperaban 6 columnas, hay {len(f)}")
        fid, iid, pid, mid, sex, pheno = f[:6]
        if fid not in fids:
            fids.append(fid)
        specs.append(IndividualSpec(
            sample=iid,
            sex={"1": "M", "2": "F"}.get(sex, "U"),
            affected=(pheno == "2"),
            phenotype_known=pheno in ("1", "2"),
            father=None if pid in ("0", "-9", ".") else pid,
            mother=None if mid in ("0", "-9", ".") else mid,
        ))
    return fids, specs


def load_ped(path: Path) -> tuple[list[str], list[IndividualSpec]]:
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"PED no encontrado: {p}")
    return parse_ped_text(p.read_text())


def write_ped(individuals: list[IndividualSpec], family_id: str = "FAM1") -> str:
    """Serializa a PED. Fenotipo desconocido se escribe como 0."""
    rows = []
    for i in individuals:
        sex = {"M": "1", "F": "2"}.get(i.sex, "0")
        pheno = ("2" if i.affected else "1") if i.phenotype_known else "0"
        rows.append("\t".join([family_id, i.sample, i.father or "0",
                               i.mother or "0", sex, pheno]))
    return "\n".join(rows) + "\n"


def validate_pedigree(
    individuals: list[IndividualSpec],
    family_ids: list[str] | None = None,
    sequenced: set[str] | None = None,
) -> PedValidation:
    """Comprobaciones estructurales. `sequenced` = muestras con BAM."""
    v = PedValidation()
    if not individuals:
        v.issues.append("pedigrí vacío")
        return v
    ids = [i.sample for i in individuals]
    dup = sorted({x for x in ids if ids.count(x) > 1})
    if dup:
        v.issues.append(f"individuos duplicados: {dup}")
    by = {i.sample: i for i in individuals}

    for i in individuals:
        for role, pid, want in (("padre", i.father, "M"), ("madre", i.mother, "F")):
            if pid is None:
                continue
            if pid == i.sample:
                v.issues.append(f"{i.sample} figura como su propio {role}")
                continue
            par = by.get(pid)
            if par is None:
                v.issues.append(
                    f"{role} de {i.sample} ('{pid}') no tiene fila en el PED")
                continue
            if par.sex == "U":
                v.warnings.append(f"sexo desconocido para {pid}, {role} de {i.sample}")
            elif par.sex != want:
                v.issues.append(
                    f"{pid} es {role} de {i.sample} pero su sexo es {par.sex}")
        if i.father and i.father == i.mother:
            v.issues.append(f"{i.sample}: padre y madre son el mismo individuo")
        if (i.father is None) != (i.mother is None):
            v.warnings.append(f"{i.sample} tiene un solo progenitor declarado")
        if not i.phenotype_known:
            v.warnings.append(f"fenotipo desconocido para {i.sample}")

    # Ciclos: nadie puede ser su propio antepasado.
    state: dict[str, int] = {}

    def visit(x: str, path: list[str]) -> None:
        if state.get(x) == 2 or x not in by:
            return
        if state.get(x) == 1:
            v.issues.append(f"ciclo en el pedigrí: {' -> '.join(path + [x])}")
            return
        state[x] = 1
        for pid in (by[x].father, by[x].mother):
            if pid:
                visit(pid, path + [x])
        state[x] = 2

    for x in ids:
        visit(x, [])

    if family_ids and len(family_ids) > 1:
        v.warnings.append(f"varias familias en el PED: {family_ids}")
    if not any(i.affected for i in individuals):
        v.warnings.append("ningún individuo afectado")
    if sequenced is not None:
        extra = sorted(sequenced - set(ids))
        if extra:
            v.issues.append(f"muestras con BAM que no están en el PED: {extra}")
        missing = sorted(i.sample for i in individuals
                         if i.affected and i.sample not in sequenced)
        if missing:
            v.warnings.append(f"afectados sin BAM: {missing}")
    return v


# ---------------------------------------------------------------------------
# Parentesco esperado (coeficientes de Wright / kinship)
# ---------------------------------------------------------------------------

def kinship_function(individuals: list[IndividualSpec]):
    """Devuelve phi(a, b): coeficiente de kinship esperado según el PED.

    phi(a, a) = (1 + F_a) / 2,  F_a = phi(padre_a, madre_a)
    phi(a, b) = (phi(padre_a, b) + phi(madre_a, b)) / 2, con a el más
    "joven" en orden topológico (no antepasado de b). Progenitor
    ausente = fundador (phi = 0). Los bucles de consanguinidad se
    recogen solos porque el cálculo recorre el grafo completo.
    Relatedness esperado = 2 * phi (0.5 padre-hijo, hermanos; más con
    consanguinidad).
    """
    by = {i.sample: i for i in individuals}
    depth: dict[str, int] = {}

    def d(x: str | None) -> int:
        if x is None or x not in by:
            return -1
        if x not in depth:
            depth[x] = 0  # protección frente a ciclos (validación los rechaza)
            i = by[x]
            depth[x] = 1 + max(d(i.father), d(i.mother))
        return depth[x]

    @lru_cache(maxsize=None)
    def phi(a: str | None, b: str | None) -> float:
        if a is None or b is None or a not in by or b not in by:
            return 0.0
        if a == b:
            i = by[a]
            return 0.5 * (1.0 + phi(i.father, i.mother))
        if d(a) < d(b):
            a, b = b, a
        i = by[a]
        return 0.5 * (phi(i.father, b) + phi(i.mother, b))

    return phi


def inbreeding(individuals: list[IndividualSpec], sample: str) -> float:
    """F esperado de un individuo = kinship de sus padres (1/16 en hijos de primos)."""
    by = {i.sample: i for i in individuals}
    i = by[sample]
    return kinship_function(individuals)(i.father, i.mother)


def relationship(individuals: list[IndividualSpec], a: str, b: str) -> str:
    """Etiqueta del parentesco declarado entre a y b."""
    by = {i.sample: i for i in individuals}
    if a == b:
        return "mismo individuo"
    ia, ib = by[a], by[b]
    if a in (ib.father, ib.mother) or b in (ia.father, ia.mother):
        return "padre/madre-hijo"
    pa = {p for p in (ia.father, ia.mother) if p}
    pb = {p for p in (ib.father, ib.mother) if p}
    shared = pa & pb
    if len(shared) == 2:
        return "hermanos completos"
    if len(shared) == 1:
        return "medio hermanos"
    return "otro parentesco" if kinship_function(individuals)(a, b) > 0 else "no emparentados"
