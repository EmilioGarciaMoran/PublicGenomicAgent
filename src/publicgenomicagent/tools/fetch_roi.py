from __future__ import annotations

from pathlib import Path

from ..env.runtime import ToolRuntime
from .base import FetchROIInput, FetchROIOutput


def fetch_roi(runtime: ToolRuntime, inp: FetchROIInput) -> FetchROIOutput:
    """Corta un BAM por región genómica (chr:start-end) usando samtools.

    Requiere que el BAM de entrada esté indexado (.bai presente).
    Genera un sub-BAM indexado en `inp.output_bam`.
    """
    in_bam = Path(inp.bam_path).expanduser().resolve()
    out_bam = Path(inp.output_bam).expanduser().resolve()
    out_bai = out_bam.with_suffix(out_bam.suffix + ".bai")

    if not in_bam.exists():
        raise FileNotFoundError(f"BAM de entrada no existe: {in_bam}")

    out_bam.parent.mkdir(parents=True, exist_ok=True)

    # samtools view -b -h <in> <region> -o <out>
    runtime.run("samtools", [
        "view", "-b", "-h",
        str(in_bam), inp.region,
        "-o", str(out_bam),
    ])

    # samtools index <out>
    runtime.run("samtools", ["index", str(out_bam)])

    bytes_written = out_bam.stat().st_size

    return FetchROIOutput(
        ok=True,
        tool="fetch_roi",
        message=f"Extraídas lecturas de {inp.region}",
        outputs={"bam": str(out_bam), "bai": str(out_bai)},
        output_bam=out_bam,
        output_bai=out_bai,
        bytes_written=bytes_written,
    )
