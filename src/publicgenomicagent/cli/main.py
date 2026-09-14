from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from publicgenomicagent.env.bootstrap import bootstrap_env
from publicgenomicagent.env.micromamba import env_exists, find_micromamba
from publicgenomicagent.env.paths import ENVS_DIR, PGA_ROOT
from publicgenomicagent.env.registry import load_registry
from publicgenomicagent.env.verify import verify_all

app = typer.Typer(help="PublicGenomicAgent — agente genómico conversacional")
env_app = typer.Typer(help="Gestión de entornos virtualizados")
app.add_typer(env_app, name="env")

console = Console()


def _registry_path() -> Path:
    # En dev: el repo está en la cwd; en prod: empaquetado en el wheel
    local = Path.cwd() / "envs" / "registry.yaml"
    if local.exists():
        return local
    raise FileNotFoundError("envs/registry.yaml no encontrado en la cwd")


@env_app.command("list")
def env_list() -> None:
    reg = load_registry(_registry_path())
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
    reg = load_registry(_registry_path())
    console.print(f"[dim]micromamba: {find_micromamba()}[/dim]")
    msg = bootstrap_env(reg, name, force=force)
    console.print(msg)


@env_app.command("verify")
def env_verify() -> None:
    reg = load_registry(_registry_path())
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
