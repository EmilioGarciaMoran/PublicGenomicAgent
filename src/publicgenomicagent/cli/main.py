from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from publicgenomicagent.env.bootstrap import bootstrap_env
from publicgenomicagent.env.micromamba import env_exists, find_micromamba
from publicgenomicagent.env.paths import ENVS_DIR, PGA_ROOT, registry_path
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.runtime import ToolRuntime
from publicgenomicagent.env.verify import verify_all

app = typer.Typer(help="PublicGenomicAgent — agente genómico conversacional")
env_app = typer.Typer(help="Gestión de entornos virtualizados")
tool_app = typer.Typer(help="Ejecutar herramientas bioinformáticas")
app.add_typer(env_app, name="env")
app.add_typer(tool_app, name="tool")

console = Console()


def _registry():
    return load_registry(registry_path())


def _runtime():
    return ToolRuntime(_registry())


# ----------------------------- env -------------------------------------

@env_app.command("list")
def env_list() -> None:
    reg = _registry()
    table = Table(title="Entornos declarados")
    table.add_column("Entorno")
    table.add_column("Estado")
    table.add_column("Manifiesto")
    table.add_column("Descripción")
    for name, spec in reg.envs.items():
        state = "[green]creado[/green]" if env_exists(name) else "[yellow]pendiente[/yellow]"
        table.add_row(name, state, str(spec.manifest), spec.description)
    console.print(table)
    console.print(f"[dim]PGA_ROOT={PGA_ROOT}[/dim]")


@env_app.command("bootstrap")
def env_bootstrap(
    name: str = typer.Argument(..., help="Nombre del entorno a crear"),
    force: bool = typer.Option(False, "--force", "-f", help="Reconstruir aunque exista"),
) -> None:
    reg = _registry()
    console.print(f"[dim]micromamba: {find_micromamba()}[/dim]")
    msg = bootstrap_env(reg, name, force=force)
    console.print(msg)


@env_app.command("verify")
def env_verify() -> None:
    reg = _registry()
    results = verify_all(reg)
    table = Table(title="Verificación de entornos y herramientas")
    table.add_column("Entorno")
    table.add_column("Herramienta")
    table.add_column("Estado")
    table.add_column("Detalle")
    for r in results:
        color = {"ok": "green", "stale_manifest": "yellow"}.get(r.status, "red")
        table.add_row(r.env, r.tool, f"[{color}]{r.status}[/{color}]", r.detail)
    console.print(table)


# ----------------------------- tool ------------------------------------

@tool_app.command("fetch-roi")
def tool_fetch_roi(
    bam: str = typer.Option(..., "--bam", "-b", help="BAM de entrada (indexado)"),
    region: str = typer.Option(..., "--region", "-r", help="Región chr:start-end"),
    out: str = typer.Option(..., "--out", "-o", help="BAM de salida"),
) -> None:
    from publicgenomicagent.tools.base import FetchROIInput
    from publicgenomicagent.tools.fetch_roi import fetch_roi

    inp = FetchROIInput(bam_path=Path(bam), region=region, output_bam=Path(out))
    result = fetch_roi(_runtime(), inp)
    console.print(
        f"[green]OK[/green] {result.output_bam} "
        f"({result.bytes_written} bytes, índice: {result.output_bai.name})"
    )


@tool_app.command("qc")
def tool_qc(
    bam: str = typer.Option(..., "--bam", "-b", help="BAM a inspeccionar"),
    level: str = typer.Option("structural", "--level", "-l",
                              help="structural | counts | deep"),
    expected_ref: str = typer.Option(None, "--expected-ref",
                                     help=".fai de referencia esperada (opcional)"),
    expected_sample: str = typer.Option(None, "--expected-sample",
                                        help="SM esperado en el header (opcional)"),
) -> None:
    from publicgenomicagent.tools.base import QCBamInput, QCLevel
    from publicgenomicagent.tools.qc import qc_bam

    try:
        qc_level = QCLevel(level)
    except ValueError:
        console.print(f"[red]Nivel inválido:[/red] {level}. Usa structural|counts|deep")
        raise typer.Exit(code=2)

    inp = QCBamInput(
        bam_path=Path(bam),
        level=qc_level,
        expected_reference_fai=Path(expected_ref) if expected_ref else None,
        expected_sample=expected_sample,
    )
    result = qc_bam(_runtime(), inp)

    color = "green" if result.passed else "red"
    console.print(f"[{color}]{result.message}[/{color}]")

    h = result.header
    console.print(f"[bold]Header[/bold]  SO={h.sort_order}  "
                  f"refs={len(h.references)}  RG={len(h.read_groups)}  SM={h.samples}")

    if result.issues:
        console.print("\n[red]Issues[/red]")
        for i in result.issues:
            console.print(f"  [red]x[/red] {i}")
    if result.warnings:
        console.print("\n[yellow]Warnings[/yellow]")
        for w in result.warnings:
            console.print(f"  [yellow]![/yellow] {w}")

    if result.counts:
        console.print("\n[bold]Counts[/bold]")
        console.print_json(data=result.counts)

    raise typer.Exit(code=0 if result.passed else 1)


@tool_app.command("call-variants")
def tool_call_variants(
    bam: str = typer.Option(..., "--bam", "-b", help="Sub-BAM de entrada"),
    reference: str = typer.Option(..., "--reference", "-f", help="FASTA de referencia"),
    out: str = typer.Option(..., "--out", "-o", help="VCF.gz de salida"),
    region: str = typer.Option(None, "--region", "-r", help="Región chr:start-end (opcional)"),
    min_qual: int = typer.Option(20, "--min-qual", help="QUAL mínimo"),
    min_dp: int = typer.Option(5, "--min-dp", help="DP mínimo por muestra"),
) -> None:
    from publicgenomicagent.tools.base import CallVariantsInput
    from publicgenomicagent.tools.call_variants import call_variants

    inp = CallVariantsInput(
        bam_path=Path(bam),
        reference_fasta=Path(reference),
        output_vcf=Path(out),
        region=region,
        min_qual=min_qual,
        min_dp=min_dp,
    )
    result = call_variants(_runtime(), inp)
    console.print(
        f"[green]OK[/green] {result.output_vcf} "
        f"({result.variants_passing} variantes passing de "
        f"{result.variants_total}, índice: {result.output_tbi.name})"
    )


@tool_app.command("compare-vcfs")
def tool_compare_vcfs(
    baseline: str = typer.Option(..., "--baseline", "-a",
                                 help="VCF baseline (p.ej. GRCh38)"),
    candidate: str = typer.Option(..., "--candidate", "-c",
                                  help="VCF candidato (p.ej. enriquecido)"),
    out: str = typer.Option(..., "--out", "-o", help="delta.vcf.gz"),
    report: str = typer.Option(None, "--report",
                               help="delta_report.json (opcional)"),
    tsv: str = typer.Option(None, "--tsv", help="delta.tsv (opcional)"),
    sample: str = typer.Option(None, "--sample", "-s",
                               help="sample name (por defecto, el primero)"),
    ground_truth: str = typer.Option(None, "--ground-truth",
                                     help="VCF ground truth (opcional)"),
) -> None:
    from publicgenomicagent.tools.base import CompareVCFsInput
    from publicgenomicagent.tools.compare_vcfs import compare_vcfs

    inp = CompareVCFsInput(
        baseline_vcf=Path(baseline),
        candidate_vcf=Path(candidate),
        output_delta_vcf=Path(out),
        output_report=Path(report) if report else None,
        output_tsv=Path(tsv) if tsv else None,
        sample=sample,
        ground_truth_vcf=Path(ground_truth) if ground_truth else None,
    )
    result = compare_vcfs(inp)
    console.print(f"[green]OK[/green] {result.delta_vcf}")
    for k, v in result.counts.items():
        console.print(f"  {k}: {v}")
    if result.metrics.get("delta_sensitivity") is not None:
        console.print(
            f"[bold]delta sensitivity:[/bold] "
            f"{result.metrics['delta_sensitivity']:+.4f}   "
            f"[bold]delta precision:[/bold] "
            f"{result.metrics['delta_precision']:+.4f}"
        )


@tool_app.command("local-pangenome")
def tool_local_pangenome(
    cohort: str = typer.Option(..., "--cohort", "-c",
                               help="VCF de cohorte con INFO/AF"),
    reference: str = typer.Option(..., "--reference", "-f",
                                  help="FASTA de referencia del ROI (+ .fai)"),
    region: str = typer.Option(..., "--region", "-r",
                               help="Región chr:start-end"),
    out: str = typer.Option(..., "--out", "-o",
                            help="Directorio de salida"),
    min_af: float = typer.Option(0.01, "--min-af",
                                 help="Umbral de frecuencia alélica"),
    af_field: str = typer.Option("AF", "--af-field",
                                 help="Campo INFO con la frecuencia"),
    contig: str = typer.Option(None, "--contig",
                               help="Nombre del contig en el FASTA (override)"),
) -> None:
    from publicgenomicagent.tools.base import LocalPangenomeInput
    from publicgenomicagent.tools.local_pangenome import local_pangenome

    inp = LocalPangenomeInput(
        cohort_vcf=Path(cohort),
        reference_fasta=Path(reference),
        region=region,
        output_dir=Path(out),
        min_af=min_af,
        af_info_field=af_field,
        contig_name=contig,
    )
    result = local_pangenome(inp)
    console.print(f"[green]OK[/green] {result.message}")
    console.print(f"  enriched: {result.enriched_fasta}")
    console.print(f"  diff VCF: {result.diff_vcf}")
    console.print(f"  report:   {result.report_json}")
    for k, v in result.counts.items():
        console.print(f"  {k}: {v}")


@tool_app.command("build-local-graph")
def tool_build_local_graph(
    reference: str = typer.Option(..., "--reference", "-f",
                                  help="FASTA del ROI (+ .fai)"),
    cohort: str = typer.Option(..., "--cohort", "-c",
                               help="cohort VCF (con AF, SVs, etc.)"),
    out: str = typer.Option(..., "--out", "-o",
                            help="Fichero .vg de salida"),
    no_alt: bool = typer.Option(False, "--no-alt-paths",
                                help="No usar -a (desactiva alt paths)"),
    max_node: int = typer.Option(32, "--max-node-length",
                                 help="Longitud máxima de nodo (-m)"),
) -> None:
    from publicgenomicagent.tools.base import BuildLocalGraphInput
    from publicgenomicagent.tools.local_graph import build_local_graph

    inp = BuildLocalGraphInput(
        reference_fasta=Path(reference),
        cohort_vcf=Path(cohort),
        output_vg=Path(out),
        include_alt_paths=not no_alt,
        max_node_length=max_node,
    )
    result = build_local_graph(_runtime(), inp)
    console.print(f"[green]OK[/green] {result.message}")
    console.print(f"  .vg: {result.graph_vg}")
    console.print(f"  .xg: {result.graph_xg}")


@tool_app.command("align-to-graph")
def tool_align_to_graph(
    graph: str = typer.Option(..., "--graph", "-g",
                              help="Fichero .vg del grafo"),
    reads: str = typer.Option(..., "--reads", "-r",
                              help="Lecturas en FASTQ"),
    out: str = typer.Option(..., "--out", "-o",
                            help="GAM de salida"),
    pack: str = typer.Option(None, "--pack",
                             help="Pack de cobertura (opcional)"),
    min_quality: int = typer.Option(5, "--min-quality",
                                    help="Calidad mínima para pack"),
) -> None:
    from publicgenomicagent.tools.base import AlignToGraphInput
    from publicgenomicagent.tools.local_graph import align_to_graph

    inp = AlignToGraphInput(
        graph_vg=Path(graph),
        reads_fastq=Path(reads),
        output_gam=Path(out),
        output_pack=Path(pack) if pack else None,
        pack_min_quality=min_quality,
    )
    result = align_to_graph(_runtime(), inp)
    console.print(f"[green]OK[/green] {result.message}")
    console.print(f"  GAM: {result.gam}")
    if result.pack:
        console.print(f"  pack: {result.pack}")


@tool_app.command("call-from-graph")
def tool_call_from_graph(
    graph: str = typer.Option(..., "--graph", "-g",
                              help="Fichero .vg del grafo"),
    xg: str = typer.Option(..., "--xg", help="Índice .xg"),
    pack: str = typer.Option(..., "--pack", help="Pack de cobertura"),
    out: str = typer.Option(..., "--out", "-o",
                            help="VCF de salida"),
) -> None:
    from publicgenomicagent.tools.base import CallFromGraphInput
    from publicgenomicagent.tools.local_graph import call_from_graph

    inp = CallFromGraphInput(
        graph_vg=Path(graph),
        graph_xg=Path(xg),
        pack=Path(pack),
        output_vcf=Path(out),
    )
    result = call_from_graph(_runtime(), inp)
    console.print(f"[green]OK[/green] {result.message}")
    console.print(f"  VCF: {result.vcf}")
    if result.vcf_tbi.name:
        console.print(f"  TBI: {result.vcf_tbi}")
