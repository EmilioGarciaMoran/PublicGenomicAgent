from __future__ import annotations

from enum import Enum
from pathlib import Path

from pydantic import BaseModel


class ToolOutput(BaseModel):
    """Resultado común de cualquier tool del sistema."""

    ok: bool
    tool: str
    message: str = ""
    outputs: dict[str, str] = {}


# --- fetch_roi ---------------------------------------------------------

class FetchROIInput(BaseModel):
    bam_path: Path
    region: str
    output_bam: Path


class FetchROIOutput(ToolOutput):
    output_bam: Path
    output_bai: Path
    bytes_written: int


# --- qc_bam ------------------------------------------------------------

class QCLevel(str, Enum):
    STRUCTURAL = "structural"  # siempre: quickcheck + índice + header
    COUNTS = "counts"          # + idxstats + flagstat
    DEEP = "deep"              # + samtools stats (aún no implementado)


class QCBamInput(BaseModel):
    bam_path: Path
    level: QCLevel = QCLevel.STRUCTURAL
    expected_reference_fai: Path | None = None
    expected_sample: str | None = None


class QCHeader(BaseModel):
    sort_order: str | None = None          # SO: coordinate | queryname | unknown
    references: list[str] = []             # @SQ SN
    read_groups: list[str] = []            # @RG ID
    samples: list[str] = []                # @RG SM (únicos)


class QCBamOutput(ToolOutput):
    passed: bool
    issues: list[str] = []
    warnings: list[str] = []
    header: QCHeader
    counts: dict[str, object] | None = None
    stats: dict[str, object] | None = None
