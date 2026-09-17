from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path

from .base import CompareVCFsInput, CompareVCFsOutput


# ---------- lectura de VCF ----------

@dataclass(frozen=True)
class VariantKey:
    chrom: str
    pos: int
    ref: str
    alt: str

    @classmethod
    def from_fields(cls, fields: list[str]) -> "VariantKey":
        # Normalización mínima para comparación: solo REF/ALT en mayúsculas.
        # CHROM se deja como viene para preservar la representación original.
        return cls(
            chrom=fields[0],
            pos=int(fields[1]),
            ref=fields[3].upper(),
            alt=fields[4].upper(),
        )


@dataclass
class VariantRecord:
    key: VariantKey
    raw_fields: list[str]
    gt: str | None      # genotipo del sample de interés, o None
    samples: list[str]


def _open_vcf(path: Path):
    """Devuelve un file handle de texto, soportando .vcf y .vcf.gz."""
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p)


def read_vcf(
    path: Path,
    sample: str | None = None,
) -> tuple[dict[VariantKey, VariantRecord], list[str], str | None]:
    """Lee un VCF y devuelve {key: VariantRecord}, [header_lines], sample_name.

    Si sample es None, usa el primer sample del header.
    """
    variants: dict[VariantKey, VariantRecord] = {}
    header: list[str] = []
    sample_name: str | None = None
    samples: list[str] = []

    with _open_vcf(path) as f:
        for line in f:
            if line.startswith("##"):
                header.append(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    raise ValueError(f"VCF sin columnas de muestra: {path}")
                samples = fields[9:]
                if sample is None:
                    sample_name = samples[0]
                else:
                    if sample not in samples:
                        raise ValueError(
                            f"sample '{sample}' no está en {path}. "
                            f"Disponibles: {samples}"
                        )
                    sample_name = sample
                header.append(line.rstrip("\n"))
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            key = VariantKey.from_fields(fields)
            gt = None
            if sample_name is not None:
                sample_idx = samples.index(sample_name)
                fmt = fields[8].split(":")
                if "GT" in fmt:
                    gt_idx = fmt.index("GT")
                    sample_field = fields[9 + sample_idx].split(":")
                    if gt_idx < len(sample_field):
                        gt = sample_field[gt_idx]
            variants[key] = VariantRecord(
                key=key, raw_fields=fields, gt=gt, samples=samples,
            )

    return variants, header, sample_name


# ---------- clasificación ----------

def classify(
    baseline: dict[VariantKey, VariantRecord],
    candidate: dict[VariantKey, VariantRecord],
) -> dict[VariantKey, str]:
    """Clasifica cada variante de la unión."""
    status: dict[VariantKey, str] = {}
    all_keys = set(baseline) | set(candidate)

    for k in all_keys:
        in_b = k in baseline
        in_c = k in candidate
        if in_b and not in_c:
            status[k] = "lost"
        elif in_c and not in_b:
            status[k] = "recovered"
        else:
            gt_b = baseline[k].gt
            gt_c = candidate[k].gt
            if gt_b == gt_c:
                status[k] = "consistent"
            else:
                status[k] = "discordant"

    return status


# ---------- escritura ----------

def write_delta_vcf(
    out_path: Path,
    header: list[str],
    sample: str,
    baseline: dict[VariantKey, VariantRecord],
    candidate: dict[VariantKey, VariantRecord],
    status: dict[VariantKey, str],
) -> None:
    """Escribe VCF con la unión de variantes anotada con PGA_STATUS."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(out_path, "wt") as f:
        # Header original (sin #CHROM)
        for line in header:
            if line.startswith("##"):
                f.write(line + "\n")
        f.write('##INFO=<ID=PGA_STATUS,Number=1,Type=String,'
                'Description="recovered|lost|consistent|discordant">\n')
        f.write('##INFO=<ID=PGA_BASELINE_GT,Number=1,Type=String,'
                'Description="GT in baseline VCF">\n')
        f.write('##INFO=<ID=PGA_CANDIDATE_GT,Number=1,Type=String,'
                'Description="GT in candidate VCF">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                f"{sample}\n")

        # Unión, ordenada por (chrom, pos)
        for k in sorted(status, key=lambda x: (x.chrom, x.pos)):
            rec = baseline.get(k) or candidate.get(k)
            assert rec is not None
            gt_b = baseline[k].gt if k in baseline else "."
            gt_c = candidate[k].gt if k in candidate else "."
            extra = f"PGA_STATUS={status[k]};PGA_BASELINE_GT={gt_b};PGA_CANDIDATE_GT={gt_c}"
            chrom, pos, vid, ref, alt = rec.raw_fields[:5]
            qual = rec.raw_fields[5] if len(rec.raw_fields) > 5 else "."
            filt = rec.raw_fields[6] if len(rec.raw_fields) > 6 else "."
            info = rec.raw_fields[7] if len(rec.raw_fields) > 7 else "."
            merged_info = f"{info};{extra}" if info and info != "." else extra
            f.write(
                f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t"
                f"{merged_info}\tGT\t{gt_c}\n"
            )


def write_report(
    out_path: Path,
    status: dict[VariantKey, str],
    baseline: dict[VariantKey, VariantRecord],
    candidate: dict[VariantKey, VariantRecord],
    baseline_sample: str | None,
    candidate_sample: str | None,
    ground_truth: dict[VariantKey, VariantRecord] | None = None,
) -> dict:
    """Contadores y métricas. Devuelve el dict reportado."""
    counts: dict[str, int] = {}
    for s in status.values():
        counts[s] = counts.get(s, 0) + 1

    report: dict[str, object] = {
        "counts": counts,
        "total_variants_in_union": len(status),
        "baseline_total": len(baseline),
        "candidate_total": len(candidate),
        "baseline_sample": baseline_sample,
        "candidate_sample": candidate_sample,
    }

    if ground_truth is not None:
        # Verdaderos positivos = variantes en ground truth detectadas por el pipeline
        # Simplificación: se compara por clave, sin genotipo por ahora
        tp_b = sum(1 for k in ground_truth if k in baseline)
        fn_b = sum(1 for k in ground_truth if k not in baseline)
        fp_b = sum(1 for k in baseline if k not in ground_truth)

        tp_c = sum(1 for k in ground_truth if k in candidate)
        fn_c = sum(1 for k in ground_truth if k not in candidate)
        fp_c = sum(1 for k in candidate if k not in ground_truth)

        def safe_div(n, d):
            return round(n / d, 4) if d > 0 else 0.0

        report["ground_truth_total"] = len(ground_truth)
        report["baseline_metrics"] = {
            "tp": tp_b, "fp": fp_b, "fn": fn_b,
            "sensitivity": safe_div(tp_b, tp_b + fn_b),
            "precision": safe_div(tp_b, tp_b + fp_b),
        }
        report["candidate_metrics"] = {
            "tp": tp_c, "fp": fp_c, "fn": fn_c,
            "sensitivity": safe_div(tp_c, tp_c + fn_c),
            "precision": safe_div(tp_c, tp_c + fp_c),
        }
        report["delta_sensitivity"] = round(
            report["candidate_metrics"]["sensitivity"]
            - report["baseline_metrics"]["sensitivity"],
            4,
        )
        report["delta_precision"] = round(
            report["candidate_metrics"]["precision"]
            - report["baseline_metrics"]["precision"],
            4,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2))
    return report


def write_tsv(
    out_path: Path,
    status: dict[VariantKey, str],
    baseline: dict[VariantKey, VariantRecord],
    candidate: dict[VariantKey, VariantRecord],
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write("chrom\tpos\tref\talt\tstatus\tgt_baseline\tgt_candidate\n")
        for k in sorted(status, key=lambda x: (x.chrom, x.pos)):
            gt_b = baseline[k].gt if k in baseline else "."
            gt_c = candidate[k].gt if k in candidate else "."
            f.write(
                f"{k.chrom}\t{k.pos}\t{k.ref}\t{k.alt}\t{status[k]}\t"
                f"{gt_b}\t{gt_c}\n"
            )


# ---------- entry point ----------

def compare_vcfs(inp: CompareVCFsInput) -> CompareVCFsOutput:
    """Compara dos VCFs y clasifica las variantes de la unión."""
    base_path = Path(inp.baseline_vcf).expanduser().resolve()
    cand_path = Path(inp.candidate_vcf).expanduser().resolve()
    delta_path = Path(inp.output_delta_vcf).expanduser().resolve()

    if not base_path.exists():
        raise FileNotFoundError(f"baseline VCF no existe: {base_path}")
    if not cand_path.exists():
        raise FileNotFoundError(f"candidate VCF no existe: {cand_path}")

    baseline, header_b, sample_b = read_vcf(base_path, sample=inp.sample)
    candidate, header_c, sample_c = read_vcf(cand_path, sample=inp.sample)

    status = classify(baseline, candidate)

    write_delta_vcf(delta_path, header_b, sample_c or sample_b or "SAMPLE",
                    baseline, candidate, status)

    report_path = None
    report = None
    if inp.output_report is not None:
        report_path = Path(inp.output_report).expanduser().resolve()
        gt = None
        if inp.ground_truth_vcf is not None:
            gt_path = Path(inp.ground_truth_vcf).expanduser().resolve()
            if gt_path.exists():
                gt, _, _ = read_vcf(gt_path)
        report = write_report(report_path, status, baseline, candidate,
                              sample_b, sample_c, ground_truth=gt)

    tsv_path = None
    if inp.output_tsv is not None:
        tsv_path = Path(inp.output_tsv).expanduser().resolve()
        write_tsv(tsv_path, status, baseline, candidate)

    counts = {}
    if report is not None:
        counts = dict(report.get("counts", {}))

    return CompareVCFsOutput(
        ok=True,
        tool="compare_vcfs",
        message=f"Comparados {len(baseline)} vs {len(candidate)}: {counts}",
        outputs={"delta_vcf": str(delta_path)},
        delta_vcf=delta_path,
        report_json=report_path,
        tsv=tsv_path,
        counts=counts,
        metrics={k: v for k, v in (report or {}).items()
                 if k not in ("counts",)},
    )
