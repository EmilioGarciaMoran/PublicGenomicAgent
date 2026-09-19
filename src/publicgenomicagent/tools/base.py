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


# --- call_variants -----------------------------------------------------

class CallVariantsInput(BaseModel):
    bam_path: Path                          # sub-BAM ya cortado por fetch_roi
    reference_fasta: Path                   # hg38.fa o referencia sintética
    output_vcf: Path                        # <out>.vcf.gz
    region: str | None = None               # opcional; si None, todo el BAM
    samples: list[str] = []                 # opcional; multi-muestra (trío)
    min_qual: int = 20
    min_dp: int = 5
    ploidy: int = 2


class CallVariantsOutput(ToolOutput):
    output_vcf: Path
    output_tbi: Path
    variants_total: int
    variants_passing: int


# --- compare_vcfs ------------------------------------------------------

class CompareVCFsInput(BaseModel):
    baseline_vcf: Path
    candidate_vcf: Path
    output_delta_vcf: Path
    output_report: Path | None = None
    output_tsv: Path | None = None
    sample: str | None = None       # si None, toma el primer sample de cada VCF
    ground_truth_vcf: Path | None = None  # opcional


class CompareVCFsOutput(ToolOutput):
    delta_vcf: Path
    report_json: Path | None = None
    tsv: Path | None = None
    counts: dict[str, int] = {}
    metrics: dict[str, object] = {}


# --- local_pangenome ---------------------------------------------------

class LocalPangenomeInput(BaseModel):
    cohort_vcf: Path
    reference_fasta: Path
    region: str                     # "chr:start-end"
    output_dir: Path
    min_af: float = 0.01
    af_info_field: str = "AF"       # campo INFO con la frecuencia
    contig_name: str | None = None  # si el FASTA usa otro nombre de contig


class LocalPangenomeOutput(ToolOutput):
    enriched_fasta: Path
    enriched_fai: Path
    diff_vcf: Path
    diff_tbi: Path
    report_json: Path
    counts: dict[str, int] = {}
