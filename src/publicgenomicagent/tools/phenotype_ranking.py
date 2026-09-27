from __future__ import annotations

import csv
import json
import shutil
import subprocess
from pathlib import Path

from .base import PhenotypeRankingInput, PhenotypeRankingOutput


# ---------------------------------------------------------------------
# Configuración
# ---------------------------------------------------------------------

DEFAULT_LIRICAL_DIR = Path.home() / ".pga" / "cache" / "lirical"


def _resolve_lirical(inp: PhenotypeRankingInput) -> tuple[Path, Path, Path]:
    """Devuelve (java_jar, data_dir, lirical_root)."""
    root = Path(inp.lirical_dir or DEFAULT_LIRICAL_DIR).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(
            f"LIRICAL no encontrado en {root}. "
            "Descarga el ZIP de https://github.com/TheJacksonLaboratory/LIRICAL/releases "
            f"y descomprímelo en {root}."
        )

    # Buscar el JAR principal (lirical-cli-*.jar)
    jars = list(root.rglob("lirical-cli-*.jar"))
    jars = [j for j in jars if "distribution" not in j.name]
    if not jars:
        raise FileNotFoundError(
            f"No se encontró lirical-cli-*.jar en {root}"
        )
    jar = jars[0]

    # El data dir vive junto al JAR
    data_dir = jar.parent / "data"
    if not data_dir.exists():
        raise FileNotFoundError(
            f"Directorio de datos de LIRICAL no encontrado: {data_dir}. "
            "Ejecuta: java -jar <jar> download -d <data_dir>"
        )

    return jar, data_dir, root


def _validate_java(java_bin: str) -> None:
    if shutil.which(java_bin) is None:
        raise FileNotFoundError(
            f"Java no encontrado: {java_bin}. "
            "LIRICAL requiere Java 17+. Instálalo con: "
            "sudo apt install openjdk-17-jre"
        )


def _build_command(
    inp: PhenotypeRankingInput,
    jar: Path,
    data_dir: Path,
    prefix: str,
) -> list[str]:
    """Construye el comando java para LIRICAL prioritize."""
    cmd = [
        inp.java_bin, "-jar", str(jar), "prioritize",
        "-p", ",".join(inp.hpo_ids),
        "-d", str(data_dir),
        "-x", prefix,
        "-o", str(inp.output_dir),
        "-f", "tsv",
    ]
    if inp.negated_hpo_ids:
        cmd += ["-n", ",".join(inp.negated_hpo_ids)]
    if inp.sex:
        cmd += ["--sex", inp.sex]
    if inp.age:
        cmd += ["--age", str(inp.age)]
    if inp.vcf:
        cmd += ["--vcf", str(Path(inp.vcf).expanduser().resolve())]
        cmd += ["--assembly", inp.assembly]
    return cmd


def _parse_tsv(tsv_path: Path, top_n: int) -> tuple[list[dict], int]:
    """Parsea el TSV de LIRICAL. Devuelve (top_candidates, total).

    El TSV tiene:
      - Líneas de comentario que empiezan por '!'
      - Una cabecera con: rank, diseaseName, diseaseCurie, pretestprob,
        posttestprob, compositeLR
      - Filas ordenadas por rank ascendente
    """
    candidates: list[dict] = []
    total = 0

    with tsv_path.open() as f:
        reader = csv.DictReader(
            (line for line in f if not line.startswith("!")),
            delimiter="\t",
        )
        for row in reader:
            total += 1
            if len(candidates) < top_n:
                candidates.append({
                    "rank": int(row["rank"]) if row.get("rank") else None,
                    "disease_name": row.get("diseaseName", ""),
                    "disease_curie": row.get("diseaseCurie", ""),
                    "pretest_prob": row.get("pretestprob", ""),
                    "posttest_prob": row.get("posttestprob", ""),
                    "composite_lr": row.get("compositeLR", ""),
                })

    return candidates, total


# ---------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------

def phenotype_ranking(inp: PhenotypeRankingInput) -> PhenotypeRankingOutput:
    """Prioriza enfermedades candidatas a partir de HPO usando LIRICAL.

    Modo phenotype-only (sin VCF): devuelve lista de enfermedades
    rankeadas por probabilidad post-test.

    Modo genotipo-aware (con VCF): añade información de variantes al
    ranking. Requiere `--assembly` (hg19 o hg38).
    """
    if not inp.hpo_ids:
        raise ValueError("hpo_ids no puede estar vacío")

    outdir = Path(inp.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    jar, data_dir, _ = _resolve_lirical(inp)
    _validate_java(inp.java_bin)

    # Prefijo único basado en un hash corto del input
    import hashlib
    key = "|".join(sorted(inp.hpo_ids)) + (str(inp.vcf) if inp.vcf else "")
    prefix = "lirical_" + hashlib.sha256(key.encode()).hexdigest()[:8]

    cmd = _build_command(inp, jar, data_dir, prefix)

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=inp.timeout_seconds,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"LIRICAL falló (exit {result.returncode}).\n"
            f"stderr:\n{result.stderr[-3000:]}"
        )

    # Localizar el TSV de salida
    tsv_path = outdir / f"{prefix}.tsv"
    if not tsv_path.exists():
        # fallback por si el prefijo cambió
        candidates = list(outdir.glob("*.tsv"))
        if not candidates:
            raise RuntimeError(
                f"LIRICAL no generó TSV en {outdir}.\n"
                f"stdout (últimas 20 líneas):\n"
                + "\n".join(result.stdout.splitlines()[-20:])
            )
        tsv_path = candidates[0]

    top_candidates, total = _parse_tsv(tsv_path, inp.top_n)

    # Report
    report_path = outdir / "phenotype_ranking_report.json"
    report_path.write_text(json.dumps({
        "hpo_ids": inp.hpo_ids,
        "negated_hpo_ids": inp.negated_hpo_ids,
        "vcf": str(inp.vcf) if inp.vcf else None,
        "assembly": inp.assembly if inp.vcf else None,
        "sex": inp.sex,
        "age": inp.age,
        "mode": "genotype_aware" if inp.vcf else "phenotype_only",
        "total_diseases_ranked": total,
        "top_n_included": len(top_candidates),
    }, indent=2))

    counts = {
        "total_diseases": total,
        "top_n": len(top_candidates),
    }

    return PhenotypeRankingOutput(
        ok=True,
        tool="phenotype_ranking",
        message=(
            f"LIRICAL: {total} enfermedades rankeadas, "
            f"top-{len(top_candidates)} incluido en el report"
        ),
        outputs={
            "ranking_tsv": str(tsv_path),
            "report": str(report_path),
        },
        ranking_tsv=tsv_path,
        top_candidates=top_candidates,
        counts=counts,
        report_json=report_path,
    )
