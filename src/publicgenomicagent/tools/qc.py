from __future__ import annotations

import shutil
from pathlib import Path

from ..env.runtime import ToolRuntime
from .base import QCHeader, QCBamInput, QCBamOutput, QCLevel


def _parse_header(runtime: ToolRuntime, bam: Path) -> QCHeader:
    """Ejecuta `samtools view -H <bam>` y parsea SO, @SQ, @RG/SM."""
    raw = runtime.capture("samtools", ["view", "-H", str(bam)])
    sort_order: str | None = None
    references: list[str] = []
    read_groups: list[str] = []
    samples: list[str] = []

    for line in raw.splitlines():
        if not line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        tag = fields[0]
        if tag == "@HD":
            for f in fields[1:]:
                if f.startswith("SO:"):
                    sort_order = f.split(":", 1)[1]
        elif tag == "@SQ":
            for f in fields[1:]:
                if f.startswith("SN:"):
                    references.append(f.split(":", 1)[1])
        elif tag == "@RG":
            for f in fields[1:]:
                if f.startswith("ID:"):
                    read_groups.append(f.split(":", 1)[1])
                elif f.startswith("SM:"):
                    sm = f.split(":", 1)[1]
                    if sm not in samples:
                        samples.append(sm)

    return QCHeader(
        sort_order=sort_order,
        references=references,
        read_groups=read_groups,
        samples=samples,
    )


def _check_index(bam: Path) -> tuple[bool, str]:
    """Devuelve (existe, ruta) del índice .bai o .csi."""
    for suffix in (".bai", ".csi"):
        cand = bam.with_suffix(bam.suffix + suffix)
        if cand.exists():
            return True, str(cand)
    # también aceptar <bam>.bai sin doble sufijo, por si acaso
    for cand in (Path(str(bam) + ".bai"), Path(str(bam) + ".csi")):
        if cand.exists():
            return True, str(cand)
    return False, ""


def _quickcheck_passed(runtime: ToolRuntime, bam: Path) -> bool:
    """samtools quickcheck devuelve exit 0 si OK, !=0 si hay problemas."""
    import subprocess
    from ..env.micromamba import bin_path

    samtools = bin_path("pga-hts", "samtools")
    if not samtools.exists():
        raise FileNotFoundError(f"samtools no encontrado en pga-hts: {samtools}")
    result = subprocess.run(
        [str(samtools), "quickcheck", "-v", str(bam)],
        capture_output=True,
        text=True,
    )
    return result.returncode == 0


def _load_reference_names(fai: Path) -> list[str]:
    """Devuelve los nombres de contig declarados en un .fai."""
    names: list[str] = []
    with fai.open() as f:
        for line in f:
            if not line.strip():
                continue
            names.append(line.split("\t", 1)[0])
    return names


def _flagstat(runtime: ToolRuntime, bam: Path) -> dict[str, object]:
    raw = runtime.capture("samtools", ["flagstat", str(bam)])
    out: dict[str, object] = {}
    for line in raw.splitlines():
        # Formato: "1234 + 0 in total (QC-passed reads + QC-failed reads)"
        parts = line.split(" ", 3)
        if len(parts) < 3:
            continue
        try:
            qc_pass = int(parts[0])
        except ValueError:
            continue
        # etiqueta aproximada
        label = " ".join(parts[3:]).split("(", 1)[0].strip() if len(parts) > 3 else line
        out[label] = qc_pass
    return out


def _idxstats(runtime: ToolRuntime, bam: Path) -> dict[str, int]:
    raw = runtime.capture("samtools", ["idxstats", str(bam)])
    out: dict[str, int] = {}
    for line in raw.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        chrom = fields[0]
        try:
            mapped = int(fields[2])
        except ValueError:
            continue
        out[chrom] = mapped
    return out


def qc_bam(runtime: ToolRuntime, inp: QCBamInput) -> QCBamOutput:
    """Control de calidad de un BAM. Estructural siempre; COUNTS opcional."""
    bam = Path(inp.bam_path).expanduser().resolve()

    issues: list[str] = []
    warnings: list[str] = []

    if not bam.exists():
        return QCBamOutput(
            ok=False,
            tool="qc_bam",
            message=f"BAM no existe: {bam}",
            passed=False,
            issues=[f"BAM no existe: {bam}"],
            header=QCHeader(),
        )

    # 1. quickcheck — integridad del fichero
    if not _quickcheck_passed(runtime, bam):
        issues.append("samtools quickcheck: fichero BAM truncado o corrupto")

    # 2. índice
    has_index, index_path = _check_index(bam)
    if not has_index:
        issues.append("no se encontró índice .bai/.csi (ejecuta: samtools index)")

    # 3. header
    header = _parse_header(runtime, bam)

    if header.sort_order not in ("coordinate", None):
        warnings.append(f"orden del BAM es '{header.sort_order}' (recomendado: coordinate)")

    if not header.read_groups:
        warnings.append("no hay @RG en el header (necesario para calling multi-muestra)")

    if header.read_groups and not header.samples:
        issues.append("@RG presente pero sin SM (sample) — calling de trío imposible")

    if len(header.samples) > 1:
        warnings.append(
            f"múltiples SM en el mismo BAM: {header.samples} "
            "(esperado 1 BAM = 1 muestra)"
        )

    # 4. referencia esperada
    if inp.expected_reference_fai is not None:
        fai = Path(inp.expected_reference_fai).expanduser().resolve()
        if not fai.exists():
            issues.append(f".fai de referencia no encontrado: {fai}")
        else:
            ref_names = _load_reference_names(fai)
            missing = [c for c in header.references if c not in ref_names]
            extra = [c for c in ref_names if c not in header.references]
            if missing:
                issues.append(
                    f"contigs del BAM ausentes en referencia: {missing[:5]}"
                    + (" ..." if len(missing) > 5 else "")
                )
            if extra and len(extra) <= 5:
                warnings.append(f"referencia tiene contigs no usados por el BAM: {extra}")

    # 5. muestra esperada
    if inp.expected_sample is not None:
        if inp.expected_sample not in header.samples:
            issues.append(
                f"SM esperado '{inp.expected_sample}' no está en header {header.samples}"
            )

    # 6. counts opcional
    counts: dict[str, object] | None = None
    if inp.level in (QCLevel.COUNTS, QCLevel.DEEP):
        if not has_index:
            warnings.append("no se puede hacer idxstats/flagstat sin índice")
        else:
            counts = {
                "flagstat": _flagstat(runtime, bam),
                "idxstats": _idxstats(runtime, bam),
            }

    if inp.level == QCLevel.DEEP:
        warnings.append("nivel DEEP no implementado aún (samtools stats pendiente)")

    passed = len(issues) == 0
    msg = "QC OK" if passed else f"QC falla con {len(issues)} problema(s)"

    return QCBamOutput(
        ok=True,
        tool="qc_bam",
        message=msg,
        outputs={"bam": str(bam), "index": index_path} if has_index else {"bam": str(bam)},
        passed=passed,
        issues=issues,
        warnings=warnings,
        header=header,
        counts=counts,
    )
