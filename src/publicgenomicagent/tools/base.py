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


# --- joint_call -------------------------------------------------------

class JointCallInput(BaseModel):
    """Variant calling conjunto (multi-sample) sobre N BAMs.

    Cada elemento de `bams` es (sample_name, bam_path). El sample_name
    debe coincidir con el campo SM del header del BAM, porque
    bcftools mpileup usa ese campo para nombrar las columnas del VCF.
    """
    bams: list[tuple[str, Path]]
    reference_fasta: Path
    output_vcf: Path
    region: str | None = None
    min_qual: int = 20
    ploidy: int = 2


class JointCallOutput(ToolOutput):
    output_vcf: Path
    output_tbi: Path
    variants_total: int
    variants_passing: int
    samples: list[str] = []


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


# --- build_local_graph ------------------------------------------------

class BuildLocalGraphInput(BaseModel):
    reference_fasta: Path
    cohort_vcf: Path
    output_vg: Path
    include_alt_paths: bool = True    # -a en vg construct
    max_node_length: int = 32         # -m
    region: str | None = None         # opcional, no usado en esta versión


class BuildLocalGraphOutput(ToolOutput):
    graph_vg: Path
    graph_xg: Path
    nodes: int = 0
    edges: int = 0


# --- align_to_graph ---------------------------------------------------

class AlignToGraphInput(BaseModel):
    graph_vg: Path
    reads_fastq: Path
    output_gam: Path
    output_pack: Path | None = None
    pack_min_quality: int = 5


class AlignToGraphOutput(ToolOutput):
    gam: Path
    pack: Path | None = None
    reads_aligned: int = 0
    reads_total: int = 0


# --- call_from_graph --------------------------------------------------

class CallFromGraphInput(BaseModel):
    graph_vg: Path
    graph_xg: Path
    pack: Path
    output_vcf: Path
    sample_name: str | None = None    # si se da, renombra SAMPLE en el VCF


class CallFromGraphOutput(ToolOutput):
    vcf: Path
    vcf_tbi: Path
    variants_total: int = 0


# --- mendelian_filter ------------------------------------------------

class IndividualSpec(BaseModel):
    """Un individuo del pedigrí, con su nombre en el VCF y sus relaciones.

    El LLM puede extraer esta información del relato clínico, o puede
    venir de un .fam, o de un JSON estructurado. La tool solo consume
    la estructura resultante.
    """
    sample: str                        # nombre en el VCF
    sex: str = "U"                     # "M" | "F" | "U"
    affected: bool = False
    father: str | None = None          # sample del padre, si aplica
    mother: str | None = None          # sample de la madre, si aplica


class MendelianFilterInput(BaseModel):
    trio_vcf: Path
    output_dir: Path
    pedigree: list[IndividualSpec]
    proband: str                       # sample del probando a analizar
    min_dp: int = 10
    min_qual: int = 20
    filter_by_affected: bool = True    # si False, dominante no exige
                                       # progenitor afectado (penetrancia
                                       # incompleta, variantes de riesgo)


class MendelianFilterOutput(ToolOutput):
    de_novo_vcf: Path
    auto_rec_hom_vcf: Path
    auto_dom_vcf: Path
    x_linked_rec_vcf: Path | None = None
    report_json: Path
    counts: dict[str, int] = {}


# --- plink_validate -------------------------------------------------

class PlinkValidateInput(BaseModel):
    trio_vcf: Path
    output_dir: Path
    pedigree: list[IndividualSpec]
    proband: str


class PlinkValidateOutput(ToolOutput):
    bed_prefix: Path
    mendel_errors: Path
    ibd_report: Path
    report_json: Path
    n_mendel_errors: int = 0


# --- CaseManifest ---------------------------------------------------

class GenomicRange(BaseModel):
    chrom: str
    start: int
    end: int
    label: str = ""    # por ejemplo, nombre del gen


class CaseManifest(BaseModel):
    """Contenedor flexible de un caso.

    No todos los casos tienen todos los componentes. El sistema debe
    degradar funcionalidad según lo que falte:

    - Sin pedigree: no se puede aplicar mendelian_filter.
    - Sin hpo_terms: no se puede priorizar por fenotipo.
    - Sin candidate_rois: el sistema debe derivarlos (por HPO o por
      búsqueda bibliográfica).

    Al menos uno de {pedigree, hpo_terms, candidate_rois} debe estar
    presente.
    """
    case_id: str
    source: str = "clinical"      # clinical | bibliography | cohort | population

    # Estructura familiar (opcional)
    pedigree: list[IndividualSpec] = []
    proband: str | None = None
    consanguinity: bool = False

    # Fenotipo (opcional)
    hpo_terms: dict[str, list[str]] = {}   # sample -> [HPO IDs]
    phenotype_text: dict[str, str] = {}    # sample -> relato libre
    affected_samples: list[str] = []

    # Regiones de interés (opcional)
    candidate_rois: list[GenomicRange] = []
    candidate_genes: list[str] = []

    # Trazabilidad
    confidence: float = 1.0
    extractor: str = ""            # "llm:model" | "manual" | "fam-file"
    confirmed_by: str | None = None

    def has_pedigree(self) -> bool:
        return len(self.pedigree) > 0 and self.proband is not None

    def has_phenotype(self) -> bool:
        return len(self.hpo_terms) > 0 or len(self.phenotype_text) > 0

    def has_rois(self) -> bool:
        return len(self.candidate_rois) > 0

    def validate_minimum(self) -> None:
        """Verifica que el manifiesto tiene al menos un componente."""
        if not (self.has_pedigree() or self.has_phenotype() or self.has_rois()):
            raise ValueError(
                "CaseManifest requiere al menos uno de: pedigree, "
                "hpo_terms/phenotype_text, o candidate_rois"
            )

    def get_parents(self, sample: str) -> tuple[str | None, str | None]:
        """Devuelve (father, mother) para un sample."""
        for ind in self.pedigree:
            if ind.sample == sample:
                return ind.father, ind.mother
        return None, None


# --- extract_hpo (RDMA) ----------------------------------------------

class ExtractHPOInput(BaseModel):
    text: str
    output_dir: Path
    backend: str = "local"                    # local | openrouter | api | azure | llama_cpp
    model_type: str = "mistral_24b"
    device: str = "auto"
    rdma_cache_dir: Path | None = None        # ~/.pga/cache/rdma
    model_cache_dir: Path | None = None       # ~/.pga/cache/models
    extractor_type: str = "retrieval"         # simple | iterative | multi | retrieval
    verifier_version: str = "v4"              # v2 | v3 | v4
    negation: bool = True
    family_history: bool = True
    skip_verification: bool = False
    top_k: int = 5


class ExtractHPOOutput(ToolOutput):
    hpo_terms_json: Path
    entities_json: Path
    hpo_ids: list[str] = []
    counts: dict[str, int] = {}
    report_json: Path


# --- phenotype_ranking (LIRICAL) -------------------------------------

class PhenotypeRankingInput(BaseModel):
    hpo_ids: list[str]
    negated_hpo_ids: list[str] = []
    output_dir: Path
    vcf: Path | None = None                # opcional, modo genotipo-aware
    assembly: str = "hg38"                 # hg19 | hg38
    sex: str = "UNKNOWN"                   # MALE | FEMALE | UNKNOWN
    age: str | None = None                 # edad en formato ISO o texto libre
    top_n: int = 50                        # candidatos a incluir en el output
    lirical_dir: Path | None = None        # ~/.pga/cache/lirical
    java_bin: str = "java"
    timeout_seconds: int = 600             # 10 minutos máximo


class PhenotypeRankingOutput(ToolOutput):
    ranking_tsv: Path
    top_candidates: list[dict] = []
    counts: dict[str, int] = {}
    report_json: Path
