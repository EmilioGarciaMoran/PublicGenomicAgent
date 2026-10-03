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
    sample_name: str = typer.Option(None, "--sample-name",
                                    help="Renombrar el sample del VCF (por defecto: SAMPLE)"),
) -> None:
    from publicgenomicagent.tools.base import CallFromGraphInput
    from publicgenomicagent.tools.local_graph import call_from_graph

    inp = CallFromGraphInput(
        graph_vg=Path(graph),
        graph_xg=Path(xg),
        pack=Path(pack),
        output_vcf=Path(out),
        sample_name=sample_name,
    )
    result = call_from_graph(_runtime(), inp)
    console.print(f"[green]OK[/green] {result.message}")
    console.print(f"  VCF: {result.vcf}")
    if result.vcf_tbi.name:
        console.print(f"  TBI: {result.vcf_tbi}")


@tool_app.command("mendelian-filter")
def tool_mendelian_filter(
    vcf: str = typer.Option(..., "--vcf", "-i",
                            help="VCF del trío"),
    out: str = typer.Option(..., "--out", "-o",
                            help="Directorio de salida"),
    pedigree_fam: str = typer.Option(None, "--pedigree-fam",
                                     help=".fam con el pedigrí"),
    proband: str = typer.Option(..., "--proband",
                                help="Sample del probando"),
    affected: str = typer.Option("", "--affected",
                                 help="Samples afectados, coma-separados"),
    filter_by_affected: bool = typer.Option(
        True, "--filter-by-affected/--no-filter-by-affected",
        help="Exigir progenitor afectado en dominante. "
             "Desactivar para penetrancia incompleta.",
    ),
    min_dp: int = typer.Option(10, "--min-dp"),
    min_qual: int = typer.Option(20, "--min-qual"),
) -> None:
    from publicgenomicagent.tools.base import MendelianFilterInput
    from publicgenomicagent.tools.mendelian import mendelian_filter
    from publicgenomicagent.tools.pedigree import load_pedigree_from_fam

    if not pedigree_fam:
        console.print("[red]Falta --pedigree-fam[/red]")
        raise typer.Exit(code=2)

    pedigree = load_pedigree_from_fam(Path(pedigree_fam))

    # Si el usuario pasa --affected, sobreescribir el flag affected
    # en cada IndividualSpec correspondiente.
    if affected:
        affected_list = [s.strip() for s in affected.split(",") if s.strip()]
        for ind in pedigree:
            ind.affected = ind.sample in affected_list

    inp = MendelianFilterInput(
        trio_vcf=Path(vcf),
        output_dir=Path(out),
        pedigree=pedigree,
        proband=proband,
        min_dp=min_dp,
        min_qual=min_qual,
        filter_by_affected=filter_by_affected,
    )

    result = mendelian_filter(inp)
    console.print(f"[green]OK[/green] {result.message}")
    for k, v in result.counts.items():
        console.print(f"  {k}: {v}")


@tool_app.command("plink-validate")
def tool_plink_validate(
    vcf: str = typer.Option(..., "--vcf", "-i",
                            help="VCF del trío"),
    out: str = typer.Option(..., "--out", "-o",
                            help="Directorio de salida"),
    pedigree_fam: str = typer.Option(..., "--pedigree-fam",
                                     help=".fam con el pedigrí"),
    proband: str = typer.Option(..., "--proband",
                                help="Sample del probando"),
) -> None:
    from publicgenomicagent.tools.base import PlinkValidateInput
    from publicgenomicagent.tools.mendelian import plink_validate
    from publicgenomicagent.tools.pedigree import load_pedigree_from_fam

    pedigree = load_pedigree_from_fam(Path(pedigree_fam))
    inp = PlinkValidateInput(
        trio_vcf=Path(vcf),
        output_dir=Path(out),
        pedigree=pedigree,
        proband=proband,
    )
    result = plink_validate(inp)
    console.print(f"[green]OK[/green] {result.message}")
    console.print(f"  mendel: {result.mendel_errors}")
    console.print(f"  ibd:    {result.ibd_report}")


@tool_app.command("extract-hpo")
def tool_extract_hpo(
    text: str = typer.Option(..., "--text", "-t",
                             help="Texto clínico (en inglés)"),
    out: str = typer.Option(..., "--out", "-o",
                            help="Directorio de salida"),
    backend: str = typer.Option("local", "--backend",
                                help="local | openrouter | api | azure | llama_cpp"),
    model: str = typer.Option("mistral_24b", "--model",
                              help="Tipo de modelo (mistral_24b, llama3_70b, ...)"),
    device: str = typer.Option("auto", "--device",
                               help="cuda:0 | cpu | auto"),
    extractor: str = typer.Option("retrieval", "--extractor",
                                  help="simple | iterative | multi | retrieval"),
    skip_verification: bool = typer.Option(
        False, "--skip-verification",
        help="Omitir la fase de verificación (más rápido, menos preciso)",
    ),
    rdma_cache: str = typer.Option(None, "--rdma-cache",
                                   help="Directorio de cache RDMA (default: ~/.pga/cache/rdma)"),
) -> None:
    from publicgenomicagent.tools.base import ExtractHPOInput
    from publicgenomicagent.tools.extract_hpo import extract_hpo

    inp = ExtractHPOInput(
        text=text,
        output_dir=Path(out),
        backend=backend,
        model_type=model,
        device=device,
        extractor_type=extractor,
        skip_verification=skip_verification,
        rdma_cache_dir=Path(rdma_cache) if rdma_cache else None,
    )
    result = extract_hpo(inp)
    console.print(f"[green]OK[/green] {result.message}")
    for k, v in result.counts.items():
        console.print(f"  {k}: {v}")
    console.print(f"  → {result.hpo_terms_json}")


@tool_app.command("phenotype-ranking")
def tool_phenotype_ranking(
    hpo: str = typer.Option(..., "--hpo", "-p",
                            help="HPO IDs separados por coma (HP:0000083,HP:0004322)"),
    out: str = typer.Option(..., "--out", "-o",
                            help="Directorio de salida"),
    negated: str = typer.Option("", "--negated",
                                help="HPO negados, separados por coma"),
    vcf: str = typer.Option(None, "--vcf",
                            help="VCF del paciente (modo genotipo-aware, opcional)"),
    assembly: str = typer.Option("hg38", "--assembly",
                                 help="hg19 | hg38 (solo con VCF)"),
    sex: str = typer.Option("UNKNOWN", "--sex",
                            help="MALE | FEMALE | UNKNOWN"),
    age: str = typer.Option(None, "--age", help="Edad del probando"),
    top_n: int = typer.Option(50, "--top-n", help="Candidatos a incluir en el report"),
    lirical_dir: str = typer.Option(None, "--lirical-dir",
                                    help="Directorio de LIRICAL (default: ~/.pga/cache/lirical)"),
) -> None:
    from publicgenomicagent.tools.base import PhenotypeRankingInput
    from publicgenomicagent.tools.phenotype_ranking import phenotype_ranking

    hpo_ids = [h.strip() for h in hpo.split(",") if h.strip()]
    negated_ids = [h.strip() for h in negated.split(",") if h.strip()]

    inp = PhenotypeRankingInput(
        hpo_ids=hpo_ids,
        negated_hpo_ids=negated_ids,
        output_dir=Path(out),
        vcf=Path(vcf) if vcf else None,
        assembly=assembly,
        sex=sex,
        age=age,
        top_n=top_n,
        lirical_dir=Path(lirical_dir) if lirical_dir else None,
    )
    result = phenotype_ranking(inp)
    console.print(f"[green]OK[/green] {result.message}")
    console.print(f"  TSV: {result.ranking_tsv}")
    console.print(f"  Top 5:")
    for c in result.top_candidates[:5]:
        console.print(
            f"    {c['rank']:>3}  {c['disease_name']}  "
            f"({c['disease_curie']})  post={c['posttest_prob']}"
        )


@tool_app.command("annotate-variants")
def tool_annotate_variants(
    vcf: str = typer.Option(..., "--vcf", "-i", help="VCF de entrada"),
    annotations_vcf: str = typer.Option(
        ..., "--annotations-vcf", "-a",
        help="VCF con anotaciones (ClinVar, gnomAD, cohorte)",
    ),
    out: str = typer.Option(..., "--out", "-o", help="VCF anotado de salida"),
    columns: str = typer.Option(
        "INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT",
        "--columns",
        help="Columnas a copiar (separadas por coma)",
    ),
    no_rename_chrs: bool = typer.Option(
        False, "--no-rename-chrs",
        help="No normalizar contigs chr1 <-> 1",
    ),
    region: str = typer.Option(None, "--region", "-r", help="Anotar solo un ROI"),
) -> None:
    """Anota un VCF con columnas de un VCF de referencia."""
    from publicgenomicagent.tools.base import AnnotateVariantsInput
    from publicgenomicagent.tools.annotate import annotate_variants

    cols = [c.strip() for c in columns.split(",") if c.strip()]
    inp = AnnotateVariantsInput(
        vcf=Path(vcf),
        annotations_vcf=Path(annotations_vcf),
        output_vcf=Path(out),
        columns=cols,
        rename_chrs=not no_rename_chrs,
        region=region,
    )
    result = annotate_variants(_runtime(), inp)
    console.print(f"[green]OK[/green] {result.message}")
    console.print(f"  VCF: {result.output_vcf}")
    console.print(f"  anotadas: {result.variants_annotated} / {result.variants_total}")


# ----------------------------- phenotype -------------------------------

phenotype_app = typer.Typer(help="Bootstrap de la capa de fenotipo")
app.add_typer(phenotype_app, name="phenotype")


@phenotype_app.command("bootstrap")
def phenotype_bootstrap(
    only_rdma: bool = typer.Option(False, "--only-rdma",
                                   help="Solo descargar vector stores de RDMA"),
    only_lirical: bool = typer.Option(False, "--only-lirical",
                                      help="Solo instalar LIRICAL"),
    with_model: bool = typer.Option(False, "--with-model",
                                    help="Descargar Mistral 24B (~14 GB, requiere GPU)"),
) -> None:
    """Descarga y configura la capa de fenotipo (RDMA + LIRICAL)."""
    import subprocess
    import sys
    from publicgenomicagent.env.paths import repo_root

    script = repo_root() / "scripts" / "phenotype_bootstrap.py"
    if not script.exists():
        console.print(f"[red]Script no encontrado:[/red] {script}")
        raise typer.Exit(code=2)

    args = [sys.executable, str(script)]
    if only_rdma:
        args.append("--only-rdma")
    if only_lirical:
        args.append("--only-lirical")
    if with_model:
        args.append("--with-model")

    result = subprocess.run(args, check=False)
    raise typer.Exit(code=result.returncode)


# ----------------------------- registry --------------------------------

@tool_app.command("list")
def tool_list() -> None:
    """Lista todas las tools registradas en el agente."""
    from publicgenomicagent.tools.registry import TOOL_REGISTRY

    table = Table(title="Tools registradas")
    table.add_column("Nombre")
    table.add_column("Runtime")
    table.add_column("Tags")
    table.add_column("Descripción")
    for name, spec in sorted(TOOL_REGISTRY.items()):
        table.add_row(
            name,
            "[green]sí[/green]" if spec.needs_runtime else "[dim]no[/dim]",
            ", ".join(spec.tags),
            spec.description,
        )
    console.print(table)


@tool_app.command("describe")
def tool_describe(
    name: str = typer.Argument(..., help="Nombre de la tool"),
) -> None:
    """Muestra el esquema Pydantic de una tool."""
    import json
    from publicgenomicagent.tools.registry import describe_tool

    try:
        d = describe_tool(name)
    except KeyError as e:
        console.print(f"[red]{e}[/red]")
        raise typer.Exit(code=2)

    console.print(f"[bold]{d['name']}[/bold] — {d['description']}")
    console.print(f"needs_runtime: {d['needs_runtime']}")
    console.print(f"tags: {', '.join(d['tags'])}")
    console.print("\n[bold]Input schema:[/bold]")
    console.print_json(json.dumps(d["input_schema"]))


# ----------------------------- run -------------------------------------

@app.command("run-trio")
def run_trio(
    father: str = typer.Option(..., "--father", help="BAM del padre"),
    mother: str = typer.Option(..., "--mother", help="BAM de la madre"),
    proband: str = typer.Option(..., "--proband", help="BAM del probando"),
    reference: str = typer.Option(..., "--reference", "-f", help="FASTA de referencia"),
    region: str = typer.Option(..., "--region", "-r", help="chr:start-end"),
    case_id: str = typer.Option("trio", "--case-id"),
    label: str = typer.Option("", "--label", help="Etiqueta del ROI (p.ej. gen)"),
    max_steps: int = typer.Option(50, "--max-steps"),
    planner_kind: str = typer.Option(
        "rule", "--planner", "-p",
        help="rule (determinista) | llm (LLM local con fallback)",
    ),
    clinvar: str = typer.Option(
        None, "--clinvar",
        help="VCF de ClinVar (opcional). Si se pasa, se anota el VCF final.",
    ),
) -> None:
    """Ejecuta el loop agéntico sobre un trío.

    Con --planner rule (default), usa el RuleBasedPlanner determinista.
    Con --planner llm, usa el LLMPlanner contra el Ollama local; si el
    LLM alucina (tool desconocida, args inválidos, stop prematuro),
    cae automáticamente al RuleBasedPlanner.
    """
    from publicgenomicagent.agent.loop import AgentLoop, build_trio_case
    from publicgenomicagent.agent.state import AgentState, SessionContext
    from publicgenomicagent.tools.base import GenomicRange

    chrom, _, coords = region.partition(":")
    start_s, _, end_s = coords.partition("-")
    roi = GenomicRange(
        chrom=chrom,
        start=int(start_s),
        end=int(end_s),
        label=label,
    )

    case = build_trio_case(
        case_id=case_id,
        father_bam=Path(father),
        mother_bam=Path(mother),
        proband_bam=Path(proband),
        roi=roi,
        reference_fasta=Path(reference),
    )

    state = AgentState(case=case)
    state.bams = {
        "father": Path(father),
        "mother": Path(mother),
        "proband": Path(proband),
    }
    state.references["hg38"] = Path(reference)

    if clinvar:
        clinvar_path = Path(clinvar)
        if not clinvar_path.exists():
            console.print(f"[red]ClinVar no existe:[/red] {clinvar_path}")
            raise typer.Exit(code=2)
        state.references["clinvar"] = clinvar_path
        console.print(f"[dim]ClinVar configurado: {clinvar_path.name}[/dim]")

    ctx = SessionContext(state=state, runtime=_runtime())
    from publicgenomicagent.agent.planner import (
        RuleBasedPlanner,
        RuleFirstLLMPlanner,
    )

    planner = None
    if planner_kind == "llm":
        # El flag "llm" usa el planner híbrido RuleFirstLLMPlanner:
        # reglas primero, LLM como desambiguador. El LLMPlanner puro
        # (LLM-first con fallback) sigue disponible en el módulo para
        # tests y usos programáticos.
        from publicgenomicagent.agent.config import load_config
        from publicgenomicagent.agent.llm_factory import build_llm_client

        cfg = load_config().llm
        console.print(
            f"[dim]planner=llm (hybrid) provider={cfg.provider} "
            f"model={cfg.model} include_paths={cfg.include_paths}[/dim]"
        )
        client = build_llm_client(cfg)
        planner = RuleFirstLLMPlanner(
            client=client,
            rule=RuleBasedPlanner(),
            include_paths=cfg.include_paths,
        )
    elif planner_kind == "rule":
        planner = RuleBasedPlanner()
    else:
        console.print(
            f"[red]--planner debe ser 'rule' o 'llm', no '{planner_kind}'[/red]"
        )
        raise typer.Exit(code=2)

    loop = AgentLoop(max_steps=max_steps, planner=planner)
    results = loop.run(ctx)

    console.print(f"[bold]Loop terminado:[/bold] {len(results)} iteraciones")
    for r in results:
        if r.executed:
            console.print(f"  [green]✓[/green] {r.step_name} → {r.tool_name}")
        else:
            console.print(f"  [yellow]·[/yellow] {r.step_name}: {r.skipped_reason}")

    console.print("\n[bold]Notas:[/bold]")
    for n in state.notes:
        console.print(f"  - {n}")

    out_json = Path("results") / f"{case_id}.session.json"
    state.to_json(out_json)
    console.print(f"\n[green]Estado guardado:[/green] {out_json}")


# ----------------------------- llm -------------------------------------

llm_app = typer.Typer(help="Cliente LLM local (Ollama / llama.cpp)")
app.add_typer(llm_app, name="llm")


@llm_app.command("ping")
def llm_ping(
    model: str = typer.Option(None, "--model", help="Override del modelo"),
    endpoint: str = typer.Option(None, "--endpoint", help="Override del endpoint"),
) -> None:
    """Hace un complete trivial al LLM configurado."""
    import time
    from publicgenomicagent.agent.config import load_config
    from publicgenomicagent.agent.llm_factory import build_llm_client

    cfg = load_config().llm
    if model:
        cfg.model = model
    if endpoint:
        cfg.endpoint = endpoint

    console.print(f"[dim]provider={cfg.provider} model={cfg.model}[/dim]")
    console.print(f"[dim]endpoint={cfg.endpoint}[/dim]")

    client = build_llm_client(cfg)
    t0 = time.perf_counter()
    try:
        out = client.complete("Responde solo 'pong'.", "ping")
    except Exception as e:  # noqa: BLE001
        console.print(f"[red]ERROR[/red] {type(e).__name__}: {e}")
        raise typer.Exit(code=1)
    dt = time.perf_counter() - t0

    console.print(f"[green]OK[/green] en {dt:.2f}s")
    console.print(f"respuesta: {out[:200]}")


@app.command("plan")
def plan(
    session_json: str = typer.Argument(..., help="Estado de sesión guardado"),
    include_paths: bool = typer.Option(
        False, "--include-paths",
        help="Enviar rutas absolutas al LLM (privacidad: por defecto no)",
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run",
        help="Solo construir el prompt, no llamar al LLM",
    ),
) -> None:
    """Carga un session.json y pide al LLMPlanner la siguiente acción."""
    import json as _json
    from publicgenomicagent.agent.config import load_config
    from publicgenomicagent.agent.llm_factory import build_llm_client
    from publicgenomicagent.agent.planner import LLMPlanner
    from publicgenomicagent.agent.prompts import render_user_prompt
    from publicgenomicagent.agent.state import AgentState

    p = Path(session_json)
    if not p.exists():
        console.print(f"[red]No existe:[/red] {p}")
        raise typer.Exit(code=2)

    state = AgentState.model_validate_json(p.read_text())

    if dry_run:
        console.print(render_user_prompt(state, include_paths=include_paths))
        raise typer.Exit(code=0)

    cfg = load_config().llm
    client = build_llm_client(cfg)
    planner = LLMPlanner(client=client, include_paths=include_paths)
    action = planner.next_action(state)

    console.print("[bold]Acción propuesta:[/bold]")
    if action is None:
        console.print("  (stop)")
    else:
        console.print(f"  tool: {action.tool_name}")
        console.print(f"  args: {action.args}")
        console.print(f"  rationale: {action.rationale}")

    console.print("\n[bold]Notas del planner:[/bold]")
    for n in state.notes:
        console.print(f"  - {n}")

# ----------------------------- report ----------------------------------

@app.command("report")
def report(
    session_json: str = typer.Argument(..., help="Estado de sesión guardado"),
    out: str = typer.Option(
        "report.html", "--out", "-o",
        help="Fichero HTML de salida",
    ),
) -> None:
    """Genera un informe HTML autocontenido a partir de un session.json."""
    from publicgenomicagent.agent.report import render_session_report

    src_path = Path(session_json)
    if not src_path.exists():
        console.print(f"[red]No existe:[/red] {src_path}")
        raise typer.Exit(code=2)

    try:
        out_path = render_session_report(src_path, Path(out))
    except Exception as e:  # noqa: BLE001
        console.print(f"[red]ERROR[/red] {type(e).__name__}: {e}")
        raise typer.Exit(code=1)

    console.print(f"[green]Informe generado:[/green] {out_path}")
    console.print(f"[dim]Ábrelo con: xdg-open {out_path}[/dim]")

# ----------------------------- igv -------------------------------------

@app.command("igv")
def igv(
    session_json: str = typer.Argument(..., help="Estado de sesión guardado"),
    out: str = typer.Option(
        None, "--out", "-o",
        help="Fichero session.xml de salida (por defecto: <case_id>.igv.xml)",
    ),
) -> None:
    """Genera un session.xml de IGV y URLs de control para el estado."""
    from publicgenomicagent.agent.igv import build_session_from_state

    src_path = Path(session_json)
    if not src_path.exists():
        console.print(f"[red]No existe:[/red] {src_path}")
        raise typer.Exit(code=2)

    info = build_session_from_state(src_path)

    # Escribir el XML
    if out is None:
        case_id = src_path.stem.replace(".session", "")
        out_path = src_path.parent / f"{case_id}.igv.xml"
    else:
        out_path = Path(out)
    out_path.write_text(info["xml"], encoding="utf-8")

    console.print(f"[green]session.xml generado:[/green] {out_path}")
    console.print(f"[dim]Ábrelo con IGV Desktop: File > Open Session[/dim]")
    console.print()
    console.print(f"[bold]Locus:[/bold] {info['locus'] or '(no definido)'}")
    console.print(f"[bold]Tracks:[/bold] {len(info['tracks'])}")
    for name, path in info["tracks"]:
        console.print(f"  - {name}  [dim]{path}[/dim]")
    console.print()
    console.print("[bold]Control URLs (IGV Desktop en localhost:60151):[/bold]")
    console.print(f"  [green]goto[/green]  {info['goto_url']}")
    for t in info["track_urls"]:
        console.print(f"  [green]load[/green]  {t['url']}")
    console.print()
    console.print(
        "[dim]Nota: las URLs de control solo funcionan si IGV Desktop "
        "está abierto en el puerto 60151 (por defecto).[/dim]"
    )

