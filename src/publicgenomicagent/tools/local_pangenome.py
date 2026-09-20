from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
from pathlib import Path

from ..env.micromamba import bin_path
from .base import LocalPangenomeInput, LocalPangenomeOutput


# ---------------------------------------------------------------------
# Helpers de referencia
# ---------------------------------------------------------------------

def _read_fasta_seq(fasta: Path) -> tuple[str, str]:
    """Devuelve (contig_name, sequence) de un FASTA de un solo contig."""
    contig = None
    parts: list[str] = []
    with fasta.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if contig is None:
                    contig = line[1:].split()[0]
                continue
            parts.append(line)
    if contig is None:
        raise ValueError(f"FASTA sin cabecera: {fasta}")
    return contig, "".join(parts)


def _parse_region(region: str) -> tuple[str, int, int]:
    """'chr2:100-200' -> ('chr2', 100, 200). 1-based inclusive."""
    if ":" not in region or "-" not in region:
        raise ValueError(f"Región mal formada: {region}")
    chrom, coords = region.split(":", 1)
    start_s, end_s = coords.split("-", 1)
    return chrom, int(start_s), int(end_s)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------
# Verificación REF
# ---------------------------------------------------------------------

def _verify_ref_matches(
    vcf_path: Path,
    fasta: Path,
    region: str,
    contig_override: str | None = None,
) -> tuple[int, list[str]]:
    """Comprueba que cada REF del VCF coincide con el FASTA.

    Devuelve (n_variantes_verificadas, lista_de_errores).

    Consideraciones:
      - El FASTA se supone con coordenadas LOCALES al ROI (1-based).
      - El VCF también debe estar en coordenadas locales.
      - region define el offset si el VCF estuviera en coordenadas genómicas.
    """
    contig, region_start, region_end = _parse_region(region)
    fasta_contig, seq = _read_fasta_seq(fasta)

    # Si el usuario pasa un override de contig, se usa para el FASTA.
    # Si no, se asume que el FASTA usa el mismo nombre que el VCF.
    fasta_contig_effective = contig_override or contig
    if fasta_contig != fasta_contig_effective:
        return 0, [
            f"FASTA declara contig '{fasta_contig}' pero se esperaba "
            f"'{fasta_contig_effective}'. Usa contig_override o corrige el FASTA."
        ]

    # Determinar offset: si el VCF tiene posiciones >= region_start,
    # asumimos que está en coordenadas genómicas y aplicamos offset.
    # Si tiene posiciones pequeñas (< length FASTA), asumimos locales.
    seq_len = len(seq)
    errors: list[str] = []
    checked = 0

    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3].upper()

            if chrom != contig:
                continue

            # Heurística de offset
            if pos >= region_start and region_start > seq_len:
                local_pos = pos - region_start + 1
            else:
                local_pos = pos

            if local_pos < 1 or local_pos + len(ref) - 1 > seq_len:
                errors.append(
                    f"{chrom}:{pos} fuera del FASTA "
                    f"(local {local_pos}, len(ref)={len(ref)}, seq_len={seq_len})"
                )
                continue

            fasta_ref = seq[local_pos - 1:local_pos - 1 + len(ref)].upper()
            if fasta_ref != ref:
                errors.append(
                    f"{chrom}:{pos} REF mismatch: "
                    f"VCF='{ref}' FASTA='{fasta_ref}'"
                )
            checked += 1

    return checked, errors


# ---------------------------------------------------------------------
# Pipeline bcftools
# ---------------------------------------------------------------------

def _bcftools() -> str:
    p = bin_path("pga-hts", "bcftools")
    if not p.exists():
        raise FileNotFoundError(f"bcftools no encontrado: {p}")
    return str(p)


def _samtools() -> str:
    p = bin_path("pga-hts", "samtools")
    if not p.exists():
        raise FileNotFoundError(f"samtools no encontrado: {p}")
    return str(p)


def _slice_cohort_to_roi(
    cohort_vcf: Path, region: str, out_vcf: Path
) -> int:
    """Corta el cohort VCF a la región. Devuelve nº de variantes."""
    subprocess.run(
        [_bcftools(), "view", "-r", region, str(cohort_vcf),
         "-Oz", "-o", str(out_vcf)],
        check=True,
    )
    subprocess.run(
        [_bcftools(), "index", "-t", str(out_vcf)],
        check=True,
    )
    result = subprocess.run(
        [_bcftools(), "view", "-H", str(out_vcf)],
        check=True, capture_output=True, text=True,
    )
    return sum(1 for l in result.stdout.splitlines() if l.strip())


def _filter_by_af(
    in_vcf: Path, out_vcf: Path, min_af: float, af_field: str
) -> int:
    """Filtra por INFO/<af_field> >= min_af. Devuelve nº de variantes."""
    expr = f"INFO/{af_field}>={min_af}"
    subprocess.run(
        [_bcftools(), "view", "-i", expr, str(in_vcf),
         "-Oz", "-o", str(out_vcf)],
        check=True,
    )
    subprocess.run([_bcftools(), "index", "-t", str(out_vcf)], check=True)
    result = subprocess.run(
        [_bcftools(), "view", "-H", str(out_vcf)],
        check=True, capture_output=True, text=True,
    )
    return sum(1 for l in result.stdout.splitlines() if l.strip())


def _exclude_svs(in_vcf: Path, out_vcf: Path) -> int:
    """Excluye todo lo que no sea SNV puro.

    El consenso lineal solo puede representar SNVs sin desplazar
    coordenadas. Los indels y SVs simbólicos cambian la longitud
    de la secuencia, lo cual rompe la alineación de las lecturas
    del BAM (que están alineadas contra GRCh38).
    """
    subprocess.run(
        [_bcftools(), "view", "-v", "snps", str(in_vcf),
         "-Oz", "-o", str(out_vcf)],
        check=True,
    )
    subprocess.run([_bcftools(), "index", "-t", str(out_vcf)], check=True)
    result = subprocess.run(
        [_bcftools(), "view", "-H", str(out_vcf)],
        check=True, capture_output=True, text=True,
    )
    return sum(1 for l in result.stdout.splitlines() if l.strip())


def _force_hom_alt(in_vcf: Path, out_vcf: Path) -> None:
    """Reescribe todos los GT a 1/1 para forzar el ALT en consensus."""
    with out_vcf.open("wb") as out_f:
        subprocess.run(
            [_bcftools(), "+setGT", str(in_vcf),
             "--output-type", "z", "--", "-t", "a", "-n", "c:1/1"],
            check=True, stdout=out_f,
        )
    subprocess.run([_bcftools(), "index", "-t", str(out_vcf)], check=True)


def _consensus(
    ref_fasta: Path, forced_vcf: Path, out_fasta: Path
) -> tuple[int, list[str]]:
    """Aplica consensus. Devuelve (n_aplicadas, warnings)."""
    result = subprocess.run(
        [_bcftools(), "consensus", "-f", str(ref_fasta), str(forced_vcf)],
        check=True, capture_output=True, text=True,
    )
    out_fasta.write_text(result.stdout)

    warnings: list[str] = []
    applied = 0
    for line in result.stderr.splitlines():
        if "Applied" in line and "variants" in line:
            # Formato: "Applied N variants"
            try:
                applied = int(line.split()[1])
            except (IndexError, ValueError):
                pass
        elif line.strip():
            warnings.append(line.strip())
    return applied, warnings


def _faidx(fasta: Path) -> Path:
    """Genera .fai con samtools faidx. Devuelve la ruta del .fai."""
    subprocess.run([_samtools(), "faidx", str(fasta)], check=True)
    return Path(str(fasta) + ".fai")


# ---------------------------------------------------------------------
# Escritura de salidas
# ---------------------------------------------------------------------

def _copy_diff_vcf(src_vcf: Path, dst_vcf: Path) -> Path:
    """Copia el VCF filtrado (forzado 1/1) como diff VCF trazable.

    Se guarda el VCF con GT forzado para que quede documentado
    exactamente qué variantes contribuyeron al consenso.
    """
    subprocess.run(
        [_bcftools(), "view", str(src_vcf), "-Oz", "-o", str(dst_vcf)],
        check=True,
    )
    subprocess.run([_bcftools(), "index", "-t", str(dst_vcf)], check=True)
    return Path(str(dst_vcf) + ".tbi")


def _write_report(
    out_path: Path,
    *,
    region: str,
    min_af: float,
    af_field: str,
    cohort_vcf: Path,
    reference_fasta: Path,
    cohort_total: int,
    after_af_filter: int,
    after_sv_filter: int,
    variants_applied: int,
    warnings: list[str],
    enriched_fasta: Path,
) -> dict:
    report = {
        "region": region,
        "min_af": min_af,
        "af_info_field": af_field,
        "cohort_vcf": str(cohort_vcf),
        "cohort_vcf_sha256": _sha256(cohort_vcf),
        "reference_fasta": str(reference_fasta),
        "reference_fasta_sha256": _sha256(reference_fasta),
        "counts": {
            "cohort_total_in_region": cohort_total,
            "after_af_filter": after_af_filter,
            "after_sv_filter": after_sv_filter,
            "variants_applied": variants_applied,
        },
        "enriched_fasta": str(enriched_fasta),
        "enriched_fasta_sha256": _sha256(enriched_fasta),
        "warnings": warnings,
    }
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    return report


# ---------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------

def local_pangenome(inp: LocalPangenomeInput) -> LocalPangenomeOutput:
    """Construye un consenso enriquecido del ROI a partir de un cohort VCF.

    Fuerza los alelos ALT de las variantes comunes (AF >= min_af) para que
    la referencia resultante incluya los alelos poblacionales comunes.
    """
    cohort = Path(inp.cohort_vcf).expanduser().resolve()
    ref = Path(inp.reference_fasta).expanduser().resolve()
    outdir = Path(inp.output_dir).expanduser().resolve()

    if not cohort.exists():
        raise FileNotFoundError(f"cohort VCF no existe: {cohort}")
    if not ref.exists():
        raise FileNotFoundError(f"FASTA de referencia no existe: {ref}")

    outdir.mkdir(parents=True, exist_ok=True)

    # Rutas de salida
    enriched_fasta = outdir / "enriched_roi.fa"
    diff_vcf = outdir / "enriched_diff.vcf.gz"
    report_json = outdir / "enrichment_report.json"

    # 1. Cortar cohort al ROI
    sliced = outdir / "_step1_sliced.vcf.gz"
    cohort_total = _slice_cohort_to_roi(cohort, inp.region, sliced)

    # 2. Filtrar por AF
    filtered = outdir / "_step2_af_filtered.vcf.gz"
    after_af = _filter_by_af(sliced, filtered, inp.min_af, inp.af_info_field)

    # 3. Excluir SVs simbólicos
    no_sv = outdir / "_step3_no_sv.vcf.gz"
    after_sv = _exclude_svs(filtered, no_sv)

    if after_sv == 0:
        raise RuntimeError(
            f"Tras filtrar por AF >= {inp.min_af}, no quedan variantes "
            f"en la región {inp.region}. Ajusta el umbral o comprueba "
            f"que el cohort VCF tiene variantes en esa zona."
        )

    # 4. Verificación REF vs FASTA
    checked, ref_errors = _verify_ref_matches(
        no_sv, ref, inp.region, contig_override=inp.contig_name,
    )
    if ref_errors:
        raise RuntimeError(
            "REF mismatch entre VCF y FASTA (abortando para evitar consenso "
            "incorrecto):\n  " + "\n  ".join(ref_errors[:10])
            + (f"\n  ... y {len(ref_errors) - 10} más" if len(ref_errors) > 10 else "")
        )

    # 5. Forzar hom-alt
    forced = outdir / "_step4_forced.vcf.gz"
    _force_hom_alt(no_sv, forced)

    # 6. Consensus
    applied, warnings = _consensus(ref, forced, enriched_fasta)

    # 7. faidx
    enriched_fai = _faidx(enriched_fasta)

    # 8. Diff VCF (copia del forzado para trazabilidad)
    diff_tbi = _copy_diff_vcf(forced, diff_vcf)

    # 9. Report
    report = _write_report(
        report_json,
        region=inp.region,
        min_af=inp.min_af,
        af_field=inp.af_info_field,
        cohort_vcf=cohort,
        reference_fasta=ref,
        cohort_total=cohort_total,
        after_af_filter=after_af,
        after_sv_filter=after_sv,
        variants_applied=applied,
        warnings=warnings,
        enriched_fasta=enriched_fasta,
    )

    # 10. Limpiar intermedios (mantener solo las 3 salidas + .fai + .tbi)
    for tmp in (sliced, filtered, no_sv, forced):
        tmp.unlink(missing_ok=True)
        Path(str(tmp) + ".tbi").unlink(missing_ok=True)

    return LocalPangenomeOutput(
        ok=True,
        tool="local_pangenome",
        message=(
            f"Consenso de {inp.region}: {applied} variantes aplicadas "
            f"(de {cohort_total} en cohorte, {after_sv} tras filtros)"
        ),
        outputs={
            "enriched_fasta": str(enriched_fasta),
            "diff_vcf": str(diff_vcf),
            "report": str(report_json),
        },
        enriched_fasta=enriched_fasta,
        enriched_fai=enriched_fai,
        diff_vcf=diff_vcf,
        diff_tbi=diff_tbi,
        report_json=report_json,
        counts=report["counts"],
    )
