"""Tabla familiar para clínicos -> PED.

El clínico rellena una tabla sencilla (TSV, se abre en Excel o LibreOffice),
con palabras en vez de códigos:

    id   relation         sex     status      bam      notes
    P1   proband          male    affected    p1.bam
    FA   father           male    unaffected  fa.bam
    P2   sister           female  affected    p2.bam

y una línea de cabecera para la consanguinidad de los padres del probando:

    parents_relationship   first cousins

PGA construye el PED: crea los progenitores y antepasados no secuenciados que
hagan falta (marcados como añadidos) y, si los padres están emparentados,
los antepasados comunes que cierran el bucle. Devuelve un resumen para que el
clínico lo confirme, con el coeficiente de consanguinidad esperado del
probando. Se aceptan términos en inglés y en español.

Para relaciones fuera del vocabulario (p.ej. el hijo de un primo), las
columnas opcionales `father` y `mother` admiten el id del progenitor.
"""
from __future__ import annotations

import csv
import getpass
import hashlib
import io
import json
import unicodedata
from datetime import datetime, timezone
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

from .base import IndividualSpec
from .pedigree import inbreeding, validate_pedigree, write_ped


def _norm(s: str | None) -> str:
    s = unicodedata.normalize("NFKD", (s or "").strip().lower())
    return " ".join("".join(c for c in s if not unicodedata.combining(c)).replace("_", " ").split())


SEX = {"male": "M", "m": "M", "man": "M", "boy": "M", "varon": "M", "hombre": "M", "h": "M",
       "nino": "M", "female": "F", "f": "F", "woman": "F", "girl": "F", "mujer": "F", "nina": "F",
       "unknown": "U", "u": "U", "desconocido": "U", "": "U"}
STATUS = {"affected": (True, True), "afectado": (True, True), "afectada": (True, True),
          "yes": (True, True), "si": (True, True),
          "unaffected": (False, True), "healthy": (False, True), "sano": (False, True),
          "sana": (False, True), "no afectado": (False, True), "no afectada": (False, True),
          "no": (False, True),
          "unknown": (False, False), "desconocido": (False, False), "": (False, False)}

# relación con el probando -> (rol canónico, sexo implícito)
RELATION = {
    "proband": ("proband", None), "probando": ("proband", None), "caso indice": ("proband", None),
    "index": ("proband", None),
    "father": ("father", "M"), "padre": ("father", "M"),
    "mother": ("mother", "F"), "madre": ("mother", "F"),
    "brother": ("sibling", "M"), "hermano": ("sibling", "M"),
    "sister": ("sibling", "F"), "hermana": ("sibling", "F"), "sibling": ("sibling", None),
    "son": ("child", "M"), "hijo": ("child", "M"),
    "daughter": ("child", "F"), "hija": ("child", "F"), "child": ("child", None),
    "paternal grandfather": ("PGF", "M"), "abuelo paterno": ("PGF", "M"),
    "paternal grandmother": ("PGM", "F"), "abuela paterna": ("PGM", "F"),
    "maternal grandfather": ("MGF", "M"), "abuelo materno": ("MGF", "M"),
    "maternal grandmother": ("MGM", "F"), "abuela materna": ("MGM", "F"),
    "paternal uncle": ("p_sib", "M"), "tio paterno": ("p_sib", "M"),
    "paternal aunt": ("p_sib", "F"), "tia paterna": ("p_sib", "F"),
    "maternal uncle": ("m_sib", "M"), "tio materno": ("m_sib", "M"),
    "maternal aunt": ("m_sib", "F"), "tia materna": ("m_sib", "F"),
    "other": ("other", None), "otro": ("other", None), "otra": ("other", None),
}
PARENTS_REL = {
    "unrelated": "unrelated", "no emparentados": "unrelated", "none": "unrelated", "no": "unrelated",
    "first cousins": "first_cousins", "primos hermanos": "first_cousins",
    "primos carnales": "first_cousins", "primos": "first_cousins",
    "second cousins": "second_cousins", "primos segundos": "second_cousins",
    "uncle niece": "uncle_niece", "uncle-niece": "uncle_niece", "tio sobrina": "uncle_niece",
    "tio-sobrina": "uncle_niece",
    "unknown": "unknown", "desconocido": "unknown", "": "unknown",
}
ADDED_LABEL = {
    "father": "father (added)", "mother": "mother (added)",
    "PGF": "paternal grandfather (added)", "PGM": "paternal grandmother (added)",
    "MGF": "maternal grandfather (added)", "MGM": "maternal grandmother (added)",
    "shared_founder_m": "shared ancestor (added)", "shared_founder_f": "shared ancestor (added)",
    "pggf": "paternal great-grandfather (added)", "pggm": "paternal great-grandmother (added)",
    "mggf": "maternal great-grandfather (added)", "mggm": "maternal great-grandmother (added)",
}
EXPECTED_F = {"unrelated": Fraction(0), "first_cousins": Fraction(1, 16),
              "second_cousins": Fraction(1, 64), "uncle_niece": Fraction(1, 8)}

TEMPLATE_EN = """\
# PublicGenomicAgent: family description. One row per relative, sequenced or not.
# Use words. PGA builds the pedigree (PED) and shows it to you for confirmation.
# relation (to the proband): proband, father, mother, brother, sister, son, daughter,
#   paternal/maternal grandfather, paternal/maternal grandmother,
#   paternal/maternal uncle, paternal/maternal aunt, other (then fill father/mother with ids)
# sex: male, female, unknown      status: affected, unaffected, unknown
# bam: file of the sequenced individual (path relative to this table), empty if not sequenced
# parents_relationship: unrelated, first cousins, second cousins, uncle-niece, unknown
parents_relationship\tfirst cousins
id\trelation\tsex\tstatus\tbam\tfather\tmother\tnotes
P1\tproband\tmale\taffected\tP1.bam\t\t\t
FA\tfather\tmale\tunaffected\tFA.bam\t\t\t
MO\tmother\tfemale\tunaffected\tMO.bam\t\t\t
P2\tsister\tfemale\taffected\tP2.bam\t\t\t
S1\tbrother\tmale\tunaffected\tS1.bam\t\t\t
"""

TEMPLATE_ES = """\
# PublicGenomicAgent: descripción de la familia. Una fila por familiar, secuenciado o no.
# Escriba con palabras. PGA construye el pedigrí (PED) y se lo muestra para confirmarlo.
# relation (con el probando): probando, padre, madre, hermano, hermana, hijo, hija,
#   abuelo paterno/materno, abuela paterna/materna, tío paterno/materno, tía paterna/materna,
#   otro (entonces indique father/mother con los id de sus padres)
# sex: varón, mujer, desconocido      status: afectado, sano, desconocido
# bam: fichero del individuo secuenciado (ruta relativa a esta tabla); vacío si no lo está
# parents_relationship: no emparentados, primos hermanos, primos segundos, tío-sobrina, desconocido
parents_relationship\tprimos hermanos
id\trelation\tsex\tstatus\tbam\tfather\tmother\tnotes
P1\tprobando\tvarón\tafectado\tP1.bam\t\t\t
FA\tpadre\tvarón\tsano\tFA.bam\t\t\t
MO\tmadre\tmujer\tsano\tMO.bam\t\t\t
P2\thermana\tmujer\tafectado\tP2.bam\t\t\t
S1\thermano\tvarón\tsano\tS1.bam\t\t\t
"""


@dataclass
class FamilyTable:
    individuals: list[IndividualSpec]
    proband: str
    bams: dict[str, Path] = field(default_factory=dict)
    relation: dict[str, str] = field(default_factory=dict)      # id -> texto original
    added: list[str] = field(default_factory=list)              # creados por PGA, no secuenciados
    parents_relationship: str = "unknown"
    expected_f: Fraction | None = None
    computed_f: float | None = None
    notes: list[str] = field(default_factory=list)
    issues: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.issues

    def ped(self, family_id: str = "FAM1") -> str:
        return write_ped(self.individuals, family_id)


def _split(text: str) -> tuple[dict[str, str], list[dict[str, str]]]:
    meta, lines = {}, []
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        lines.append(raw)
    delim = "\t" if any("\t" in l for l in lines) else ("," if any("," in l for l in lines) else ";")
    rows = list(csv.reader(io.StringIO("\n".join(lines)), delimiter=delim))
    header_i = next((i for i, r in enumerate(rows) if r and _norm(r[0]) == "id"), None)
    if header_i is None:
        raise ValueError("no encuentro la fila de cabecera (debe empezar por 'id')")
    for r in rows[:header_i]:
        if len(r) >= 2:
            meta[_norm(r[0])] = r[1]
    header = [_norm(h) for h in rows[header_i]]
    data = [dict(zip(header, (c.strip() for c in r))) for r in rows[header_i + 1:] if any(c.strip() for c in r)]
    return meta, data


def parse_family_table(text: str, base_dir: Path | None = None) -> FamilyTable:
    meta, data = _split(text)
    issues, warnings, notes = [], [], []
    by: dict[str, IndividualSpec] = {}
    role: dict[str, str] = {}
    relation_txt: dict[str, str] = {}
    bams: dict[str, Path] = {}
    for n, row in enumerate(data, 1):
        iid = row.get("id", "").strip()
        if not iid:
            issues.append(f"fila {n}: falta el id")
            continue
        if iid in by:
            issues.append(f"id repetido: {iid}")
            continue
        rel_raw = row.get("relation", "")
        rel = RELATION.get(_norm(rel_raw))
        if rel is None:
            issues.append(f"{iid}: relación no reconocida '{rel_raw}' (use la lista de la plantilla u 'other')")
            continue
        sex_raw = row.get("sex", "")
        sex = SEX.get(_norm(sex_raw))
        if sex is None:
            issues.append(f"{iid}: sexo no reconocido '{sex_raw}'")
            continue
        if rel[1] and sex == "U":
            sex = rel[1]
            notes.append(f"{iid}: sexo deducido de la relación ({rel_raw})")
        elif rel[1] and sex != rel[1]:
            issues.append(f"{iid}: '{rel_raw}' no concuerda con el sexo '{sex_raw}'")
        st_raw = row.get("status", "")
        st = STATUS.get(_norm(st_raw))
        if st is None:
            issues.append(f"{iid}: estado no reconocido '{st_raw}' (affected/unaffected/unknown)")
            continue
        by[iid] = IndividualSpec(sample=iid, sex=sex, affected=st[0], phenotype_known=st[1],
                                 father=row.get("father") or None, mother=row.get("mother") or None)
        role[iid], relation_txt[iid] = rel[0], rel_raw
        if row.get("bam"):
            p = Path(row["bam"]).expanduser()
            bams[iid] = (Path(base_dir) / p if base_dir and not p.is_absolute() else p).resolve()

    probands = [i for i, r in role.items() if r == "proband"]
    if len(probands) != 1:
        issues.append(f"debe haber exactamente un probando, hay {len(probands)}")
        return FamilyTable(list(by.values()), probands[0] if probands else "", bams, relation_txt,
                           issues=issues, warnings=warnings, notes=notes)
    prob = probands[0]
    added: list[str] = []

    def one(r: str) -> str | None:
        ids = [i for i, x in role.items() if x == r]
        if len(ids) > 1:
            issues.append(f"más de un '{r}': {ids}")
        return ids[0] if ids else None

    def ensure(r: str, sex: str, label: str) -> str:
        have = one(r)
        if have:
            return have
        new = label
        while new in by:
            new += "_"
        by[new] = IndividualSpec(sample=new, sex=sex, affected=False, phenotype_known=False)
        role[new], relation_txt[new] = r, ADDED_LABEL.get(r, f"{r} (added)")
        added.append(new)
        return new

    needs_parents = any(r in ("sibling", "PGF", "PGM", "MGF", "MGM", "p_sib", "m_sib") for r in role.values())
    rel_parents = PARENTS_REL.get(_norm(meta.get("parents relationship", "")))
    if rel_parents is None:
        issues.append(f"parents_relationship no reconocido: '{meta.get('parents relationship')}'")
        rel_parents = "unknown"
    fa = one("father") or (ensure("father", "M", "FATHER") if needs_parents or rel_parents not in ("unrelated", "unknown") else None)
    mo = one("mother") or (ensure("mother", "F", "MOTHER") if needs_parents or rel_parents not in ("unrelated", "unknown") else None)
    by[prob].father, by[prob].mother = by[prob].father or fa, by[prob].mother or mo
    for i, r in role.items():
        if r == "sibling":
            by[i].father, by[i].mother = by[i].father or fa, by[i].mother or mo
        elif r == "child":
            if by[i].father is None and by[i].mother is None:
                (setattr(by[i], "father", prob) if by[prob].sex == "M" else setattr(by[i], "mother", prob))
                warnings.append(f"{i}: el otro progenitor no figura en la tabla")

    def grandparents(side: str) -> tuple[str, str]:
        child = fa if side == "P" else mo
        gf = ensure(f"{side}GF", "M", f"{side}GF")
        gm = ensure(f"{side}GM", "F", f"{side}GM")
        by[child].father, by[child].mother = by[child].father or gf, by[child].mother or gm
        return gf, gm

    for side, sib in (("P", "p_sib"), ("M", "m_sib")):
        if any(r in (f"{side}GF", f"{side}GM", sib) for r in role.values()):
            gf, gm = grandparents(side)
            for i, r in role.items():
                if r == sib:
                    by[i].father, by[i].mother = by[i].father or gf, by[i].mother or gm

    if rel_parents in ("first_cousins", "second_cousins", "uncle_niece"):
        pgf, _ = grandparents("P")
        mgf, mgm = grandparents("M")
        if rel_parents == "first_cousins":
            # abuelo paterno y abuela materna son hermanos (hijos de los fundadores comunes)
            sf, sm = ensure("shared_founder_m", "M", "FOUNDER1"), ensure("shared_founder_f", "F", "FOUNDER2")
            for x in (pgf, mgm):
                by[x].father, by[x].mother = by[x].father or sf, by[x].mother or sm
            notes.append("padres primos hermanos: se asume que el abuelo paterno y la abuela materna "
                         "son hermanos (el parentesco autosómico no depende de esta elección; para "
                         "herencia ligada al X indique los enlaces exactos con father/mother)")
        elif rel_parents == "second_cousins":
            a = ensure("pggf", "M", "PGGF")
            b = ensure("mggm", "F", "MGGM")
            by[pgf].father = by[pgf].father or a
            by[mgm].mother = by[mgm].mother or b
            for x, other_r, other_sex, lab in ((pgf, "pggm", "F", "PGGM"), (mgm, "mggf", "M", "MGGF")):
                o = ensure(other_r, other_sex, lab)
                if other_sex == "F":
                    by[x].mother = by[x].mother or o
                else:
                    by[x].father = by[x].father or o
            sf, sm = ensure("shared_founder_m", "M", "FOUNDER1"), ensure("shared_founder_f", "F", "FOUNDER2")
            for x in (a, b):
                by[x].father, by[x].mother = by[x].father or sf, by[x].mother or sm
            notes.append("padres primos segundos: se asume que un bisabuelo paterno y una bisabuela "
                         "materna son hermanos")
        else:  # tío-sobrina: el abuelo materno es hermano del padre
            pgf2, pgm2 = by[fa].father, by[fa].mother
            by[mgf].father, by[mgf].mother = by[mgf].father or pgf2, by[mgf].mother or pgm2
            notes.append("padres tío-sobrina: se asume que el padre es hermano del abuelo materno")

    inds = list(by.values())
    v = validate_pedigree(inds, sequenced=set(bams))
    issues += v.issues
    warnings += [w for w in v.warnings
                 if not any(w.endswith(f"para {a}") for a in added)]   # los añadidos no tienen fenotipo
    for i, p in bams.items():
        if not p.exists():
            warnings.append(f"{i}: BAM no encontrado ({p})")
    exp = EXPECTED_F.get(rel_parents)
    comp = None
    if not issues:
        comp = inbreeding(inds, prob)
        if exp is not None and abs(comp - float(exp)) > 1e-9:
            issues.append(f"consanguinidad calculada F={comp:.4f} no coincide con la declarada "
                          f"({rel_parents}, F={float(exp):.4f}): revise father/mother")
    if added:
        notes.append(f"añadidos por PGA (no secuenciados, fenotipo desconocido): {', '.join(added)}")
    return FamilyTable(inds, prob, bams, relation_txt, added, rel_parents, exp, comp,
                       notes, issues, warnings)


def load_family_table(path: Path) -> FamilyTable:
    p = Path(path).expanduser().resolve()
    return parse_family_table(p.read_text(encoding="utf-8-sig"), base_dir=p.parent)


def summary(ft: FamilyTable) -> str:
    """Resumen legible para que el clínico confirme antes de analizar."""
    st = lambda i: ("affected" if i.affected else "unaffected") if i.phenotype_known else "unknown"  # noqa: E731
    w = max(len(i.sample) for i in ft.individuals)
    lines = [f"{'id':<{w}}  {'relation':<28} sex  {'status':<10} sequenced  father / mother"]
    for i in ft.individuals:
        lines.append(f"{i.sample:<{w}}  {ft.relation.get(i.sample, ''):<28} {i.sex:<4} {st(i):<10} "
                     f"{'yes' if i.sample in ft.bams else 'no':<10} {i.father or '-'} / {i.mother or '-'}")
    lines.append("")
    lines.append(f"Proband: {ft.proband}. Parents: {ft.parents_relationship.replace('_', ' ')}.")
    if ft.computed_f is not None:
        frac = Fraction(ft.computed_f).limit_denominator(1024)
        lines.append(f"Expected inbreeding coefficient of the proband: F = {frac} ({ft.computed_f:.4f}).")
    for n in ft.notes:
        lines.append(f"Note: {n}")
    for wn in ft.warnings:
        lines.append(f"Warning: {wn}")
    for e in ft.issues:
        lines.append(f"ERROR: {e}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Aprobación: el análisis solo corre con una tabla aprobada y sin cambios
# ---------------------------------------------------------------------------

def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def review_path(table: Path) -> Path:
    t = Path(table)
    return t.with_name(t.stem + ".check.json")


def record_review(table: Path, drawing: bool) -> dict:
    """Registra que se generó la página de revisión de esta versión exacta de la tabla."""
    table = Path(table).expanduser().resolve()
    rec = {"table_sha256": _sha256(table.read_bytes()), "drawing": drawing,
           "reviewed_at": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    review_path(table).write_text(json.dumps(rec, indent=2))
    return rec


def approval_path(table: Path) -> Path:
    t = Path(table)
    return t.with_name(t.stem + ".approval.json")


def ped_path(table: Path) -> Path:
    t = Path(table)
    return t.with_name(t.stem + ".ped")


def approve(table: Path, by: str | None = None, family_id: str = "FAM1") -> dict:
    """Escribe el PED y el registro de aprobación (hashes de la tabla y del PED)."""
    table = Path(table).expanduser().resolve()
    rv = review_path(table)
    if not rv.exists() or json.loads(rv.read_text()).get("table_sha256") != _sha256(table.read_bytes()):
        raise ValueError(f"revise primero esta versión de la tabla: pga family check {table.name}")
    ft = load_family_table(table)
    if not ft.ok:
        raise ValueError("la tabla tiene errores: " + "; ".join(ft.issues))
    ped = ft.ped(family_id)
    ped_path(table).write_text(ped)
    rec = {"table": table.name, "table_sha256": _sha256(table.read_bytes()),
           "ped": ped_path(table).name, "ped_sha256": _sha256(ped.encode()),
           "proband": ft.proband, "parents_relationship": ft.parents_relationship,
           "expected_inbreeding_proband": ft.computed_f,
           "added_by_pga": ft.added, "reviewed_with_drawing": json.loads(rv.read_text())["drawing"],
           "approved_by": by or getpass.getuser(),
           "approved_at": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    approval_path(table).write_text(json.dumps(rec, indent=2, ensure_ascii=False))
    return rec


def check_approval(table: Path) -> dict:
    """Devuelve el registro si la tabla actual está aprobada; si no, ValueError con el motivo."""
    table = Path(table).expanduser().resolve()
    ap = approval_path(table)
    if not ap.exists():
        raise ValueError(f"tabla sin aprobar: revise con 'pga family check {table.name}' y apruebe con "
                         f"'pga family approve {table.name} --by NOMBRE'")
    rec = json.loads(ap.read_text())
    if rec.get("table_sha256") != _sha256(table.read_bytes()):
        raise ValueError("la tabla ha cambiado desde su aprobación: revísela y apruébela de nuevo")
    pp = ped_path(table)
    if not pp.exists() or rec.get("ped_sha256") != _sha256(pp.read_bytes()):
        raise ValueError("el PED aprobado falta o ha cambiado: apruebe de nuevo")
    return rec


def check_page(ft: FamilyTable, svg_text: str | None, table_name: str) -> str:
    """Página HTML autocontenida para la confirmación del clínico."""
    import html as _h
    legend = ("Red: affected · grey: phenotype unknown · green outline: sequenced · arrow: proband · "
              "double line: consanguineous couple · “(added)”: created by PGA, not sequenced")
    drawing = svg_text or ('<p class="warn">Drawing not available (install the R environment: '
                           '<code>pga env bootstrap pga-r</code>).</p>')
    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Family check</title>
<style>body{{font-family:system-ui,sans-serif;max-width:980px;margin:24px auto;padding:0 16px;color:#1d2433}}
pre{{background:#f4f6f9;padding:12px;overflow-x:auto;font-size:13px}} .warn{{color:#b42318}}
.box{{border:1px solid #d0d5dd;border-radius:8px;padding:12px;margin:16px 0}} svg{{max-width:100%;height:auto}}</style>
</head><body><h1>Family check: {_h.escape(table_name)}</h1>
<p>Please check the drawing and the table. If they are right, approve with
<code>pga family approve {_h.escape(table_name)} --by YOUR_NAME</code>. Nothing is analysed before approval.</p>
<div class="box">{drawing}<p><small>{legend}</small></p></div>
<h2>As understood by PGA</h2><pre>{_h.escape(summary(ft))}</pre>
<h2>Pedigree file that will be used (PED)</h2><pre>{_h.escape(ft.ped())}</pre>
</body></html>"""
