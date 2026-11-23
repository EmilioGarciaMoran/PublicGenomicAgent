"""Dibujo del pedigrí con Pedixplorer (R, Bioconductor; sucesor de kinship2).

Se usa para que el clínico confirme la familia antes de analizar: afectados
en rojo, fenotipo desconocido en gris, contorno verde = secuenciado, flecha =
probando, doble línea = pareja consanguínea. Los individuos añadidos por PGA
llevan "(added)" en su etiqueta.
"""
from __future__ import annotations

import csv
import subprocess
from pathlib import Path

from ..env.micromamba import bin_path
from .family_table import FamilyTable

R_SCRIPT = r"""
args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(Pedixplorer))
df <- read.delim(args[1], stringsAsFactors = FALSE, na.strings = c("NA", ""))
df$proband <- as.logical(df$proband)
ped <- Pedigree(df)
svg(args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
plot(ped, aff_mark = FALSE, label = "relation", cex = 0.85, legend = FALSE,
     title = args[5], title_cex = 1.3)
invisible(dev.off())
"""


def rscript() -> Path:
    return bin_path("pga-r", "Rscript")


def draw_pedigree(ft: FamilyTable, out_svg: Path, title: str = "Family: draft for confirmation") -> Path:
    """Escribe el SVG del pedigrí. Requiere el entorno pga-r (pga env bootstrap pga-r)."""
    r = rscript()
    if not r.exists():
        raise FileNotFoundError("entorno pga-r no instalado: pga env bootstrap pga-r")
    out_svg = Path(out_svg).resolve()
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    tsv = out_svg.with_suffix(".plot.tsv")
    script = out_svg.with_suffix(".plot.R")
    script.write_text(R_SCRIPT)
    with tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "dadid", "momid", "sex", "affection", "avail", "proband", "relation"])
        for i in ft.individuals:
            w.writerow([i.sample, i.father or "NA", i.mother or "NA",
                        {"M": 1, "F": 2}.get(i.sex, 3),
                        (1 if i.affected else 0) if i.phenotype_known else "NA",
                        1 if i.sample in ft.bams else 0,
                        "TRUE" if i.sample == ft.proband else "FALSE",
                        ft.relation.get(i.sample, "")])
    gens = 2 + sum(1 for i in ft.individuals if i.father or i.mother) // 3
    width = max(7.0, 1.2 * len(ft.individuals) ** 0.8)
    subprocess.run([str(r), str(script), str(tsv), str(out_svg), f"{width:.1f}",
                    f"{min(10.0, 2.2 * gens):.1f}", title], check=True, capture_output=True, text=True)
    tsv.unlink(); script.unlink()
    return out_svg
