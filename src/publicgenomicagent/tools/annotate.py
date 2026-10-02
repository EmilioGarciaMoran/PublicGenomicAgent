"""VCF annotation with external VCFs (ClinVar, gnomAD, cohorts).

Wraps `bcftools annotate` to copy columns from a reference VCF
(by default: ClinVar) into a query VCF, matching by
(CHROM, POS, REF, ALT).

Handles the common case where the annotation VCF uses a
different contig naming convention (Ensembl "1" vs UCSC "chr1",
which is what ClinVar ships). When `rename_chrs=True`, contigs
are normalized to the query VCF before the join.
"""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

from ..env.micromamba import bin_path
from ..env.runtime import ToolRuntime
from .base import AnnotateVariantsInput, AnnotateVariantsOutput


def _count_variants(vcf_gz: Path) -> int:
    n = 0
    with gzip.open(vcf_gz, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                n += 1
    return n


def _count_annotated(vcf_gz: Path, columns: list[str]) -> int:
    """Cuenta variantes con al menos una columna no vacía.

    Se apoya en el header ##INFO para saber qué campos
    esperar. Si una columna no aparece en INFO, se ignora.
    """
    keys = []
    for col in columns:
        # "INFO/CLNSIG" -> "CLNSIG"
        if col.startswith("INFO/"):
            keys.append(col.split("/", 1)[1])
        elif col.startswith("FORMAT/"):
            keys.append(col.split("/", 1)[1])
    if not keys:
        return 0

    annotated = 0
    with gzip.open(vcf_gz, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = fields[7]
            if any(f"{k}=" in info for k in keys):
                annotated += 1
    return annotated


def _rename_chrs_map(header_path: Path) -> dict[str, str]:
    """Construye un mapa {contig_sin_chr: contig_con_chr}.

    Se usa para normalizar el annotations_vcf al estilo del
    query VCF (chr1, chr2, ...).
    """
    mapping = {}
    with header_path.open() as f:
        for line in f:
            if line.startswith("##contig="):
                # ##contig=<ID=1,length=...>
                id_part = line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">\n")
                if not id_part.startswith("chr"):
                    mapping[id_part] = f"chr{id_part}"
                else:
                    mapping[id_part] = id_part[3:]
    return mapping


def annotate_variants(
    runtime: ToolRuntime,
    inp: AnnotateVariantsInput,
) -> AnnotateVariantsOutput:
    """Copia columnas de un VCF de anotaciones al VCF de entrada.

    Usa `bcftools annotate` con `--rename-chrs` cuando sea
    necesario, para que coincidan las nomenclaturas de contigs.
    """
    vcf = Path(inp.vcf).expanduser().resolve()
    ann = Path(inp.annotations_vcf).expanduser().resolve()
    out_vcf = Path(inp.output_vcf).expanduser().resolve()
    out_tbi = Path(str(out_vcf) + ".tbi")

    if not vcf.exists():
        raise FileNotFoundError(f"VCF de entrada no existe: {vcf}")
    if not ann.exists():
        raise FileNotFoundError(f"VCF de anotaciones no existe: {ann}")
    if out_vcf.exists() and not inp.force:
        raise FileExistsError(f"Salida ya existe: {out_vcf}")

    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    bcftools = bin_path("pga-hts", "bcftools")
    if not bcftools.exists():
        raise FileNotFoundError(f"bcftools no encontrado: {bcftools}")

    # Normalizar contigs del ANNOTATION al estilo del QUERY.
    #
    # IMPORTANTE: no aplicamos --rename-chrs al comando principal
    # porque bcftools lo aplica al QUERY (VCF de entrada), no al
    # annotations. Lo que queremos es dejar el query intacto
    # (que mantiene "chr1", "chr2", ...) y renombrar el annotation
    # ("1", "2", ...) a nuestro estilo.
    #
    # Para ello:
    #   1. Leemos el header del query.
    #   2. Leemos el header del annotation.
    #   3. Si sus contigs difieren, creamos una versión temporal
    #      del annotation con --rename-chrs y usamos esa versión.
    #   4. El comando final NO lleva --rename-chrs.
    ann_to_use = ann
    tmp_ann: Path | None = None
    chrs_file: Path | None = None

    if inp.rename_chrs:
        # Header del query
        hdr_q = subprocess.run(
            [str(bcftools), "view", "-h", str(vcf)],
            capture_output=True, text=True, check=True,
        ).stdout
        # Header del annotation
        hdr_a = subprocess.run(
            [str(bcftools), "view", "-h", str(ann)],
            capture_output=True, text=True, check=True,
        ).stdout

        def _contigs(header: str) -> set[str]:
            cs = set()
            for line in header.splitlines():
                if line.startswith("##contig="):
                    cid = line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">")
                    cs.add(cid)
            return cs

        q_contigs = _contigs(hdr_q)
        a_contigs = _contigs(hdr_a)

        # Si el annotation tiene contigs sin "chr" y el query con "chr",
        # normalizamos el annotation.
        if any(not c.startswith("chr") for c in a_contigs) and any(
            c.startswith("chr") for c in q_contigs
        ):
            chrs_file = out_vcf.with_suffix(".rename_chrs.tmp")
            tmp_ann = out_vcf.with_suffix(".normalized_ann.vcf.gz")
            with chrs_file.open("w") as f:
                for c in sorted(a_contigs):
                    if not c.startswith("chr"):
                        f.write(f"{c}\tchr{c}\n")
                    else:
                        f.write(f"{c}\t{c}\n")
            subprocess.run(
                [str(bcftools), "annotate", "--rename-chrs",
                 str(chrs_file), "-Oz", "-o", str(tmp_ann), str(ann)],
                check=True, capture_output=True,
            )
            subprocess.run(
                [str(bcftools), "index", "-t", str(tmp_ann)],
                check=True, capture_output=True,
            )
            ann_to_use = tmp_ann

    # Construir el comando final (ya con ann_to_use resuelto)
    cmd: list[str] = [str(bcftools), "annotate"]
    cmd += ["-a", str(ann_to_use)]
    cmd += ["-c", ",".join(inp.columns)]

    if inp.region:
        cmd += ["-r", inp.region]

    cmd += ["-Oz", "-o", str(out_vcf), str(vcf)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    finally:
        for tmp in (chrs_file, tmp_ann):
            if tmp and tmp.exists():
                try:
                    tmp.unlink()
                except OSError:
                    pass
        # también limpiar el índice del tmp_ann si existe
        if tmp_ann:
            tbi = Path(str(tmp_ann) + ".tbi")
            if tbi.exists():
                try:
                    tbi.unlink()
                except OSError:
                    pass

    if result.returncode != 0:
        raise RuntimeError(
            f"bcftools annotate falló (exit {result.returncode}):\n"
            f"{result.stderr}"
        )

    runtime.run("bcftools", ["index", "-t", str(out_vcf)])

    total = _count_variants(out_vcf)
    annotated = _count_annotated(out_vcf, inp.columns)

    return AnnotateVariantsOutput(
        ok=True,
        tool="annotate_variants",
        message=f"{annotated} de {total} variantes anotadas",
        outputs={"vcf": str(out_vcf), "tbi": str(out_tbi)},
        output_vcf=out_vcf,
        output_tbi=out_tbi,
        variants_total=total,
        variants_annotated=annotated,
        columns_copied=list(inp.columns),
    )

