"""Registro lógico de tools del agente.

Distinto de `envs/registry.yaml`, que mapea herramientas binarias a
entornos micromamba. Este registro mapea **tools de dominio** (fetch_roi,
qc_bam, mendelian_filter, ...) a:

  - su modelo Pydantic de input
  - su modelo Pydantic de output
  - la función Python que la implementa
  - si necesita un ToolRuntime como primer argumento

El agente usa `dispatch(tool_name, **kwargs)` para ejecutar cualquier
tool sin conocer su firma exacta. La CLI usa `list_tools()` y
`describe_tool(name)` para introspección.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

from pydantic import BaseModel

from ..env.runtime import ToolRuntime
from .base import (
    AlignToGraphInput,
    AlignToGraphOutput,
    BuildLocalGraphInput,
    BuildLocalGraphOutput,
    CallFromGraphInput,
    CallFromGraphOutput,
    CallVariantsInput,
    CallVariantsOutput,
    CompareVCFsInput,
    CompareVCFsOutput,
    ExtractHPOInput,
    ExtractHPOOutput,
    FetchROIInput,
    FetchROIOutput,
    LocalPangenomeInput,
    LocalPangenomeOutput,
    MendelianFilterInput,
    MendelianFilterOutput,
    PhenotypeRankingInput,
    PhenotypeRankingOutput,
    PlinkValidateInput,
    PlinkValidateOutput,
    QCBamInput,
    QCBamOutput,
    ToolOutput,
)
from .call_variants import call_variants
from .compare_vcfs import compare_vcfs
from .extract_hpo import extract_hpo
from .fetch_roi import fetch_roi
from .local_graph import align_to_graph, build_local_graph, call_from_graph
from .local_pangenome import local_pangenome
from .mendelian import mendelian_filter, plink_validate
from .phenotype_ranking import phenotype_ranking
from .qc import qc_bam


@dataclass(frozen=True)
class ToolSpec:
    """Metadatos de una tool registrada."""

    name: str
    func: Callable[..., ToolOutput]
    input_model: type[BaseModel]
    output_model: type[BaseModel]
    needs_runtime: bool
    description: str
    tags: tuple[str, ...] = field(default_factory=tuple)

    def input_schema(self) -> dict[str, Any]:
        """JSON Schema del input, para el prompt del LLM o la CLI."""
        return self.input_model.model_json_schema()

    def output_schema(self) -> dict[str, Any]:
        return self.output_model.model_json_schema()


TOOL_REGISTRY: dict[str, ToolSpec] = {
    # --- HTS ----------------------------------------------------------
    "fetch_roi": ToolSpec(
        name="fetch_roi",
        func=fetch_roi,
        input_model=FetchROIInput,
        output_model=FetchROIOutput,
        needs_runtime=True,
        description="Corta un BAM indexado por región genómica (chr:start-end).",
        tags=("hts", "roi"),
    ),
    "qc_bam": ToolSpec(
        name="qc_bam",
        func=qc_bam,
        input_model=QCBamInput,
        output_model=QCBamOutput,
        needs_runtime=True,
        description="QC de un BAM: header, índice, orden, muestras, counts.",
        tags=("hts", "qc"),
    ),
    "call_variants": ToolSpec(
        name="call_variants",
        func=call_variants,
        input_model=CallVariantsInput,
        output_model=CallVariantsOutput,
        needs_runtime=True,
        description="Variant calling sobre un sub-BAM (bcftools mpileup+call).",
        tags=("hts", "variants"),
    ),
    "compare_vcfs": ToolSpec(
        name="compare_vcfs",
        func=compare_vcfs,
        input_model=CompareVCFsInput,
        output_model=CompareVCFsOutput,
        needs_runtime=False,
        description="Diff entre VCF baseline y candidato, con métricas opcionales.",
        tags=("variants", "diff"),
    ),

    # --- Grafos locales / pangenoma -----------------------------------
    "local_pangenome": ToolSpec(
        name="local_pangenome",
        func=local_pangenome,
        input_model=LocalPangenomeInput,
        output_model=LocalPangenomeOutput,
        needs_runtime=False,
        description="Enriquece una referencia local con alelos de cohorte (AF).",
        tags=("pangenome", "reference"),
    ),
    "build_local_graph": ToolSpec(
        name="build_local_graph",
        func=build_local_graph,
        input_model=BuildLocalGraphInput,
        output_model=BuildLocalGraphOutput,
        needs_runtime=True,
        description="Construye un grafo vg a partir de referencia + cohorte VCF.",
        tags=("pangenome", "graph"),
    ),
    "align_to_graph": ToolSpec(
        name="align_to_graph",
        func=align_to_graph,
        input_model=AlignToGraphInput,
        output_model=AlignToGraphOutput,
        needs_runtime=True,
        description="Alinea lecturas FASTQ contra un grafo vg (vg giraffe).",
        tags=("pangenome", "graph", "align"),
    ),
    "call_from_graph": ToolSpec(
        name="call_from_graph",
        func=call_from_graph,
        input_model=CallFromGraphInput,
        output_model=CallFromGraphOutput,
        needs_runtime=True,
        description="Llama variantes desde un grafo + pack (vg call).",
        tags=("pangenome", "graph", "variants"),
    ),

    # --- Mendelian / pedigrí ------------------------------------------
    "mendelian_filter": ToolSpec(
        name="mendelian_filter",
        func=mendelian_filter,
        input_model=MendelianFilterInput,
        output_model=MendelianFilterOutput,
        needs_runtime=False,
        description="Filtros mendelianos: de novo, recesivo hom, dominante, X.",
        tags=("family", "mendelian"),
    ),
    "plink_validate": ToolSpec(
        name="plink_validate",
        func=plink_validate,
        input_model=PlinkValidateInput,
        output_model=PlinkValidateOutput,
        needs_runtime=False,
        description="Validación con PLINK: errores Mendel, IBD, sexo.",
        tags=("family", "qc"),
    ),

    # --- Fenotipo -----------------------------------------------------
    "extract_hpo": ToolSpec(
        name="extract_hpo",
        func=extract_hpo,
        input_model=ExtractHPOInput,
        output_model=ExtractHPOOutput,
        needs_runtime=False,
        description="Extrae términos HPO de texto clínico con RDMA.",
        tags=("phenotype", "nlp"),
    ),
    "phenotype_ranking": ToolSpec(
        name="phenotype_ranking",
        func=phenotype_ranking,
        input_model=PhenotypeRankingInput,
        output_model=PhenotypeRankingOutput,
        needs_runtime=False,
        description="Ranking de candidatos con LIRICAL (HPO ± VCF).",
        tags=("phenotype", "ranking"),
    ),
}


def list_tools() -> list[str]:
    """Nombres de todas las tools registradas, ordenados."""
    return sorted(TOOL_REGISTRY.keys())


def get_tool(name: str) -> ToolSpec:
    if name not in TOOL_REGISTRY:
        raise KeyError(
            f"Tool '{name}' no registrada. "
            f"Disponibles: {', '.join(list_tools())}"
        )
    return TOOL_REGISTRY[name]


def describe_tool(name: str) -> dict[str, Any]:
    """Descripción de una tool, útil para CLI / prompt del LLM."""
    spec = get_tool(name)
    return {
        "name": spec.name,
        "description": spec.description,
        "tags": list(spec.tags),
        "needs_runtime": spec.needs_runtime,
        "input_schema": spec.input_schema(),
        "output_schema": spec.output_schema(),
    }


def dispatch(
    tool_name: str,
    runtime: ToolRuntime | None = None,
    *,
    validate_only: bool = False,
    **kwargs: Any,
) -> ToolOutput | BaseModel:
    """Ejecuta una tool por nombre, con validación Pydantic automática.

    Args:
        tool_name: nombre registrado (p.ej. "fetch_roi").
        runtime: ToolRuntime, requerido si la tool lo necesita.
        validate_only: si True, solo construye y valida el input,
                       sin ejecutar la tool. Útil para tests y dry-run.
        **kwargs: campos del modelo de input.

    Returns:
        El output Pydantic de la tool, o el input validado si validate_only.

    Raises:
        KeyError: tool no registrada.
        ValidationError: input inválido según el modelo Pydantic.
        ValueError: tool requiere runtime y no se ha pasado.
    """
    spec = get_tool(tool_name)

    # Pydantic v2 valida y lanza ValidationError directamente.
    # No lo re-empaquetamos: el original ya incluye el modelo, los
    # campos faltantes y las rutas. Re-lanzarlo con `model=` rompe
    # porque ValidationError no acepta ese kwarg en su constructor.
    inp = spec.input_model(**kwargs)

    if validate_only:
        return inp

    if spec.needs_runtime:
        if runtime is None:
            raise ValueError(
                f"La tool '{tool_name}' requiere un ToolRuntime, "
                f"pero se pasó None."
            )
        return spec.func(runtime, inp)

    return spec.func(inp)
