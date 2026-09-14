from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel


class ToolOutput(BaseModel):
    """Resultado común de cualquier tool del sistema."""

    ok: bool
    tool: str
    message: str = ""
    outputs: dict[str, str] = {}


class FetchROIInput(BaseModel):
    bam_path: Path
    region: str
    output_bam: Path


class FetchROIOutput(ToolOutput):
    output_bam: Path
    output_bai: Path
    bytes_written: int
