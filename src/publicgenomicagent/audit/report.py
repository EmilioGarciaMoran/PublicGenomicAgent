"""Render de un informe de auditoría: humano (Rich) o máquina (dict)."""
from __future__ import annotations

from typing import Any

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .reader import Session
from .verify import verify_session


def to_dict(session: Session) -> dict[str, Any]:
    """Representación serializable de la sesión, para --json."""
    st = session.state
    return {
        "path": str(session.path),
        "case_id": st.case.case_id,
        "source": st.case.source,
        "proband": st.case.proband,
        "consanguinity": st.case.consanguinity,
        "tool_calls": [
            {
                "tool_name": c.tool_name,
                "input": c.input,
                "output": c.output,
                "ok": c.ok,
                "duration_ms": c.duration_ms,
                "timestamp": c.timestamp.isoformat(),
                "error": c.error,
            }
            for c in st.tool_calls
        ],
        "notes": list(st.notes),
        "outputs_keys": sorted(st.outputs.keys()),
        "n_tool_calls": len(st.tool_calls),
        "n_notes": len(st.notes),
    }


def render_human(
    session: Session,
    console: Console,
    verbose: bool = False,
    with_verify: bool = False,
) -> None:
    """Imprime el informe en la consola dada."""
    st = session.state

    header = Table.grid(padding=(0, 1))
    header.add_column(style="bold cyan", justify="right")
    header.add_column()
    header.add_row("Sesión:", str(session.path.parent.name))
    header.add_row("Case:", st.case.case_id)
    header.add_row("Source:", st.case.source)
    if st.case.proband:
        header.add_row("Proband:", st.case.proband)
    header.add_row(
        "Tools:", f"{len(st.tool_calls)} ejecutadas, "
                  f"{len([c for c in st.tool_calls if not c.ok])} fallidas"
    )
    header.add_row("Notas:", str(len(st.notes)))

    console.print(Panel(header, title="Audit", border_style="cyan"))

    if st.tool_calls:
        t = Table(title="Tool calls", show_lines=False)
        t.add_column("#", justify="right", style="dim")
        t.add_column("tool", style="bold")
        t.add_column("ok", justify="center")
        t.add_column("dur (ms)", justify="right")
        t.add_column("args", overflow="fold")

        for i, c in enumerate(st.tool_calls, 1):
            mark = "[green]✓[/green]" if c.ok else "[red]✗[/red]"
            args_repr = ", ".join(
                f"{k}={_short(v)}" for k, v in sorted(c.input.items())
            )
            t.add_row(str(i), c.tool_name, mark, str(c.duration_ms), args_repr)
            if verbose and c.error:
                t.add_row("", f"[red]{c.error}[/red]", "", "", "")
            if verbose and c.output:
                t.add_row(
                    "", f"[dim]{_short(c.output, limit=200)}[/dim]", "", "", ""
                )

        console.print(t)
    else:
        console.print("[dim]Sin tool calls registradas.[/dim]")

    if st.notes:
        console.print()
        console.print("[bold]Notas:[/bold]")
        for n in st.notes:
            console.print(f"  [dim]·[/dim] {n}")

    if with_verify:
        issues = verify_session(session)
        console.print()
        if not issues:
            console.print("[green]✓ Integridad OK[/green]")
        else:
            console.print(f"[red]✗ {len(issues)} problemas:[/red]")
            for issue in issues:
                console.print(f"  [red]·[/red] {issue}")


def _short(v: Any, limit: int = 60) -> str:
    s = str(v)
    if len(s) > limit:
        return s[: limit - 1] + "…"
    return s
