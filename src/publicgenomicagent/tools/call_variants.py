from __future__ import annotations

import subprocess
from pathlib import Path

from ..env.micromamba import bin_path
from ..env.runtime import ToolRuntime
from .base import CallVariantsInput, CallVariantsOutput


def _count_variants(vcf_gz: Path) -> tuple[int, int]:
    """Devuelve (total, passing) contando líneas no-header.

    'passing' se calcula con el filtro de QUAL/DP que ya aplicamos
    en el pipeline, así que a este nivel total == passing en la
    mayoría de casos. Se deja el segundo conteo para futuras fases.
    """
    total = 0
    passing = 0
    try:
        import gzip
        with gzip.open(vcf_gz, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                total += 1
                passing += 1
    except FileNotFoundError:
        pass
    return total, passing


def call_variants(runtime: ToolRuntime, inp: CallVariantsInput) -> CallVariantsOutput:
    """Genera un VCF filtrado a partir de un (sub-)BAM y una referencia.

    Usa: bcftools mpileup | bcftools call -mv | bcftools filter
    Todo dentro del entorno pga-hts, sin pasar por micromamba run.
    """
    bam = Path(inp.bam_path).expanduser().resolve()
    ref = Path(inp.reference_fasta).expanduser().resolve()
    out_vcf = Path(inp.output_vcf).expanduser().resolve()
    out_tbi = Path(str(out_vcf) + ".tbi")

    if not bam.exists():
        raise FileNotFoundError(f"BAM no existe: {bam}")
    if not ref.exists():
        raise FileNotFoundError(f"Referencia no existe: {ref}")

    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    bcftools = bin_path("pga-hts", "bcftools")
    if not bcftools.exists():
        raise FileNotFoundError(f"bcftools no encontrado: {bcftools}")

    # 1. bcftools mpileup | bcftools call -mv
    #    Canalizamos con subprocess.Popen.
    mpileup_cmd: list[str] = [
        str(bcftools), "mpileup",
        "-f", str(ref),
        "-a", "FORMAT/DP,FORMAT/AD",   # anotar DP y AD por muestra
    ]
    if inp.region:
        mpileup_cmd += ["-r", inp.region]
    if inp.samples:
        for s in inp.samples:
            mpileup_cmd += ["-s", s]
    mpileup_cmd.append(str(bam))

    call_cmd: list[str] = [
        str(bcftools), "call",
        "-mv",
        "-Ou",              # BCF sin comprimir para pipe
        "--ploidy", str(inp.ploidy),
    ]
    filter_cmd: list[str] = [
        str(bcftools), "filter",
        "-i", f"QUAL>={inp.min_qual} && FORMAT/DP>={inp.min_dp}",
        "-Oz",
        "-o", str(out_vcf),
    ]

    p1 = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(call_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p3 = subprocess.Popen(filter_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Cerrar el stdout del padre para que p1 reciba SIGPIPE si p2 muere
    if p1.stdout:
        p1.stdout.close()
    if p2.stdout:
        p2.stdout.close()

    _, err3 = p3.communicate()
    _, err2 = p2.communicate()
    _, err1 = p1.communicate()

    for name, code, err in (("mpileup", p1.returncode, err1),
                            ("call", p2.returncode, err2),
                            ("filter", p3.returncode, err3)):
        if code != 0:
            raise RuntimeError(
                f"bcftools {name} falló (exit {code}):\n{err.decode(errors='replace')}"
            )

    # 2. Indexar
    runtime.run("bcftools", ["index", "-t", str(out_vcf)])

    total, passing = _count_variants(out_vcf)

    return CallVariantsOutput(
        ok=True,
        tool="call_variants",
        message=f"{passing} variantes passing de {total} totales",
        outputs={"vcf": str(out_vcf), "tbi": str(out_tbi)},
        output_vcf=out_vcf,
        output_tbi=out_tbi,
        variants_total=total,
        variants_passing=passing,
    )
