"""Fenotipo × genotipo: enfermedades de LIRICAL -> genes -> candidatas familiares.

Sin red ni bases nuevas: usa los ficheros que LIRICAL ya descarga en su
directorio de datos.

- mim2gene_medgen (NCBI): filas `phenotype` enlazan un MIM de enfermedad con
  un GeneID.
- hgnc_complete_set.txt (HGNC): GeneID (entrez_id) -> símbolo oficial, y
  símbolos previos o alias -> símbolo oficial (las etiquetas de ROI pueden
  usar nombres antiguos).

`concordance()` no combina las dos evidencias en una cifra: pone junto a cada
candidata que segrega la enfermedad mejor situada por fenotipo de su gen.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class DiseaseGenes:
    omim_to_genes: dict[str, set[str]] = field(default_factory=dict)   # "OMIM:256100" -> {"NPHP1"}
    alias_to_symbol: dict[str, str] = field(default_factory=dict)      # mayúsculas -> símbolo oficial

    def symbol(self, name: str | None) -> str | None:
        if not name:
            return None
        return self.alias_to_symbol.get(name.strip().upper())


def load_disease_genes(data_dir: Path) -> DiseaseGenes:
    data_dir = Path(data_dir)
    entrez, alias = {}, {}
    with (data_dir / "hgnc_complete_set.txt").open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym = (row.get("symbol") or "").strip()
            if not sym:
                continue
            alias.setdefault(sym.upper(), sym)
            if row.get("entrez_id"):
                entrez[row["entrez_id"].strip()] = sym
            for col in ("prev_symbol", "alias_symbol"):
                for a in (row.get(col) or "").strip('"').split("|"):
                    if a.strip():
                        alias.setdefault(a.strip().upper(), sym)   # el oficial gana siempre
    omim: dict[str, set[str]] = {}
    with (data_dir / "mim2gene_medgen").open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 3 or c[2] != "phenotype" or c[1] in ("-", ""):
                continue
            sym = entrez.get(c[1])
            if sym:
                omim.setdefault(f"OMIM:{c[0]}", set()).add(sym)
    return DiseaseGenes(omim_to_genes=omim, alias_to_symbol=alias)


def read_lirical_tsv(path: Path) -> list[dict]:
    rows = []
    with Path(path).open(encoding="utf-8") as f:
        lines = [l for l in f if not l.startswith("!")]
    for row in csv.DictReader(lines, delimiter="\t"):
        try:
            rows.append({"rank": int(row["rank"]), "disease_name": row.get("diseaseName", ""),
                         "disease_curie": row.get("diseaseCurie", ""),
                         "posttest_prob": row.get("posttestprob", ""),
                         "composite_lr": row.get("compositeLR", "")})
        except (KeyError, ValueError):
            continue
    return rows


def concordance(candidates: list[dict], ranking: list[dict], dg: DiseaseGenes,
                top_n: int = 20) -> dict:
    """candidates: [{"gene", "model", "variant"}]. ranking: filas de LIRICAL (todas)."""
    lr_count: dict[str, int] = {}
    for r in ranking:
        lr_count[r["composite_lr"]] = lr_count.get(r["composite_lr"], 0) + 1
    best: dict[str, dict] = {}
    for r in sorted(ranking, key=lambda x: x["rank"]):
        for g in dg.omim_to_genes.get(r["disease_curie"], ()):
            best.setdefault(g, {**r, "tied_with": lr_count.get(r["composite_lr"], 1) - 1})
    rows, cand_genes = [], set()
    for c in candidates:
        sym = dg.symbol(c.get("gene"))
        if sym:
            cand_genes.add(sym)
        b = best.get(sym) if sym else None
        rows.append({**c, "gene_symbol": sym,
                     "gene_status": ("not_a_gene" if sym is None else
                                     "no_ranked_disease" if b is None else "ranked"),
                     "best_disease": b})
    rows.sort(key=lambda r: (r["best_disease"] or {}).get("rank", 10**9))
    top_without = []
    for r in sorted(ranking, key=lambda x: x["rank"]):
        if r["rank"] > top_n:
            break
        genes = sorted(dg.omim_to_genes.get(r["disease_curie"], ()))
        missing = [g for g in genes if g not in cand_genes]
        if missing:
            top_without.append({**r, "genes": missing})
    return {"candidates": rows, "top_diseases_without_candidate": top_without, "top_n": top_n,
            "n_ranked_diseases": len(ranking)}
