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


@tool_app.command("fetch-roi")
def tool_fetch_roi(
    bam: str = typer.Option(..., "--bam", "-b", help="BAM de entrada (indexado)"),
    region: str = typer.Option(..., "--region", "-r", help="Región chr:start-end"),
    out: str = typer.Option(..., "--out", "-o", help="BAM de salida"),
) -> None:
    from publicgenomicagent.tools.base import FetchROIInput
    from publicgenomicagent.tools.fetch_roi import fetch_roi

    reg = _registry()
    runtime = ToolRuntime(reg)
    inp = FetchROIInput(bam_path=Path(bam), region=region, output_bam=Path(out))
    result = fetch_roi(runtime, inp)
    console.print(
        f"[green]OK[/green] {result.output_bam} "
        f"({result.bytes_written} bytes, índice: {result.output_bai.name})"
    )
