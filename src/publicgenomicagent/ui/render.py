"""Traducción Event -> Rich. Sin estado, sin lógica de negocio."""
from __future__ import annotations

from rich.console import Console

from ..agent.events import Event


def render(event: Event, console: Console, verbose: bool = False) -> None:
    k = event.kind
    p = event.payload

    if k == "step_start":
        console.print(f"[dim]· {p.get('tool', '?')}[/dim]", end=" ")

    elif k == "tool_call":
        if verbose:
            console.print(
                f"[dim]→ call {p.get('name', '?')} {p.get('args', {})}[/dim]"
            )

    elif k == "tool_result":
        ok = bool(p.get("ok", True))
        mark = "[green]✓[/green]" if ok else "[red]✗[/red]"
        console.print(mark)
        if verbose and p.get("error"):
            console.print(f"[dim red]  {p['error']}[/dim red]")

    elif k == "note":
        if verbose:
            console.print(f"[dim]· {p.get('text', '')}[/dim]")

    elif k == "final":
        console.print("[dim]— fin —[/dim]")

    elif k == "error":
        console.print(f"[red]error:[/red] {p.get('message', '')}")
