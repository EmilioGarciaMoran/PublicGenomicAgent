"""Detección de regiones de homocigosidad aparente (ROH) con bcftools.

Esta tool NO es un ROH caller genómico. Los paneles de genes y los
VCFs restringidos a regiones pequeñas no tienen densidad de SNPs
suficiente para estimar F_ROH global ni mapear ROH largos (>8 Mb).

La tool ejecuta `bcftools roh` sobre el VCF de entrada y devuelve
regiones de homocigosidad aparente dentro de la zona cubierta, como
heurística de priorización de variantes recesivas homocigotas.

El output incluye un campo `caveat` obligatorio que restringe la
interpretación clínica del resultado.
"""
from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

from ..env.runtime import ToolRuntime
from .base import (
    DetectHomozygousRegionsInput,
    DetectHomozygousRegionsOutput,
    HomozygousRegion,
)


CAVEAT = (
    "Estimación restringida al VCF de entrada. No interpretable como "
    "F_ROH genómico ni como coeficiente de consanguinidad."
)


# ---------------------------------------------------------------------------
# Parseo de la salida de bcftools roh -O r
# ---------------------------------------------------------------------------

# Líneas de datos: "RG <sample> <chrom> <start> <end> <len_bp> <n_markers> <quality>"
# Comentarios: empiezan por "#"
_ROH_DATA_LINE = re.compile(r"^RG\t")


def _parse_roh_output(text: str, sample: str) -> list[HomozygousRegion]:
    """Parsea la salida de `bcftools roh -O r`.

    Ignora líneas de comentario (empiezan por '#') y la cabecera.
    Extrae cada línea 'RG' como una región. El sample viene implícito
    porque consultamos roh con `-s <sample>`, pero lo pasamos para
    asignarlo a cada región.
    """
    regions: list[HomozygousRegion] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        if not _ROH_DATA_LINE.match(line):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        try:
            chrom = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            length_bp = int(fields[5])
            n_markers = int(fields[6])
            quality = float(fields[7])
        except (ValueError, IndexError):
            continue

        regions.append(
            HomozygousRegion(
                sample=sample,
                chrom=chrom,
                start=start,
                end=end,
                length_bp=length_bp,
                length_kb=length_bp / 1000.0,
                n_markers=n_markers,
                quality=quality,
            )
        )
    return regions


# ---------------------------------------------------------------------------
# Filtros
# ---------------------------------------------------------------------------

def _filter_regions(
    regions: list[HomozygousRegion],
    min_length_kb: int,
    min_snps: int,
    min_quality: float,
) -> list[HomozygousRegion]:
    """Aplica umbrales conservadores para descartar regiones espurias."""
    min_length_bp = min_length_kb * 1000
    return [
        r
        for r in regions
        if r.length_bp >= min_length_bp
        and r.n_markers >= min_snps
        and r.quality >= min_quality
    ]


# ---------------------------------------------------------------------------
# Ejecución por sample
# ---------------------------------------------------------------------------

def _list_samples_from_vcf(vcf: Path, runtime: ToolRuntime) -> list[str]:
    """Devuelve la lista de samples del VCF usando bcftools query."""
    result = runtime.run(
        "bcftools",
        ["query", "-l", str(vcf)],
        capture_output=True,
        text=True,
    )
    return [s.strip() for s in result.stdout.splitlines() if s.strip()]


def _run_roh_for_sample(
    vcf: Path,
    sample: str,
    output_dir: Path,
    runtime: ToolRuntime,
    skip_indels: bool,
) -> str:
    """Ejecuta `bcftools roh -O r` para un sample. Devuelve la ruta del .txt."""
    out_path = output_dir / f"roh.{sample}.txt"
    args: list[str] = [
        "roh",
        "-s", sample,
        "-e", "-",          # estimar AF del propio VCF
        "-O", "r",          # salida como regiones
        "-o", str(out_path),
    ]
    if skip_indels:
        args.append("-I")
    args.append(str(vcf))

    runtime.run("bcftools", args, capture_output=True, text=True)
    return str(out_path)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def detect_homozygous_regions(
    inp: DetectHomozygousRegionsInput,
    runtime: ToolRuntime,
) -> DetectHomozygousRegionsOutput:
    """Detecta regiones de homocigosidad aparente en un VCF."""
    vcf = Path(inp.vcf).expanduser().resolve()
    if not vcf.exists():
        raise FileNotFoundError(f"VCF no encontrado: {vcf}")

    output_dir = Path(inp.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    samples: list[str]
    if inp.samples:
        samples = list(inp.samples)
    else:
        samples = _list_samples_from_vcf(vcf, runtime)

    if not samples:
        raise ValueError(f"El VCF {vcf} no contiene samples detectables")

    # Ejecutar roh por cada sample y agregar
    all_regions: list[HomozygousRegion] = []
    per_sample_counts: dict[str, int] = {}

    for sample in samples:
        out_path = _run_roh_for_sample(
            vcf=vcf,
            sample=sample,
            output_dir=output_dir,
            runtime=runtime,
            skip_indels=inp.skip_indels,
        )
        text = Path(out_path).read_text(encoding="utf-8")
        regions = _parse_roh_output(text, sample)
        regions = _filter_regions(
            regions,
            min_length_kb=inp.min_length_kb,
            min_snps=inp.min_snps,
            min_quality=inp.min_quality,
        )
        all_regions.extend(regions)
        per_sample_counts[sample] = len(regions)

    # Métricas agregadas
    total_bp = sum(r.length_bp for r in all_regions)
    total_mb = total_bp / 1_000_000.0
    # f_roh_proxy: fracción de Mb cubiertos por ROH respecto a Mb totales
    # Si el VCF es single-sample, es la del sample. Si multi, es la suma
    # de todos los samples dividida por el número de samples (proxy).
    f_roh_proxy = total_mb / max(1, len(samples))

    # Escribir resumen
    summary_path = output_dir / "roh_summary.json"
    summary = {
        "vcf": str(vcf),
        "samples_analyzed": samples,
        "per_sample_counts": per_sample_counts,
        "total_regions": len(all_regions),
        "total_mb": round(total_mb, 4),
        "f_roh_proxy": round(f_roh_proxy, 6),
        "caveat": CAVEAT,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "thresholds": {
            "min_length_kb": inp.min_length_kb,
            "min_snps": inp.min_snps,
            "min_quality": inp.min_quality,
            "skip_indels": inp.skip_indels,
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Output
    if not all_regions:
        message = (
            f"Sin regiones de homocigosidad aparente sobre "
            f"{len(samples)} sample(s) con los umbrales actuales "
            f"(min_length_kb={inp.min_length_kb}, min_snps={inp.min_snps})"
        )
    else:
        message = (
            f"{len(all_regions)} region(es) en {len(samples)} sample(s), "
            f"{total_mb:.2f} Mb totales"
        )

        # Escribir las regiones como JSON aparte; en outputs solo paths (str)
    regions_path = output_dir / "roh_regions.json"
    regions_path.write_text(
        json.dumps(
            [r.model_dump(mode="json") for r in all_regions],
            indent=2,
        ),
        encoding="utf-8",
    )

    return DetectHomozygousRegionsOutput(
        ok=True,
        tool="detect_homozygous_regions",
        message=message,
        outputs={
            "summary": str(summary_path),
            "regions": str(regions_path),
        },
        regions=all_regions,
        samples_analyzed=samples,
        total_regions=len(all_regions),
        total_mb=total_mb,
        f_roh_proxy=f_roh_proxy,
        caveat=CAVEAT,
    )
