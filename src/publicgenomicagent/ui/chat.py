"""REPL conversacional sobre AgentLoop.stream()."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from rich.console import Console

from ..agent.events import Event
from ..agent.loop import AgentLoop
from ..agent.planner import RuleBasedPlanner
from ..agent.state import AgentState, SessionContext
from ..tools.base import CaseManifest
from .render import render

console = Console()

SESSIONS_DIR = Path("sessions")


def _build_context(
    case_path: Path | None,
    use_llm: bool,
) -> tuple[SessionContext, object | None]:
    """Construye SessionContext + planner.

    Sin --case, se usa un CaseManifest mínimo ("adhoc") sin BAMs ni ROIs.
    Sin --llm, el planner es RuleBasedPlanner (determinista, sin red).
    """
    if case_path is not None:
        import json
        raw = json.loads(case_path.read_text(encoding="utf-8"))
        manifest = CaseManifest.model_validate(raw)
    else:
        manifest = CaseManifest(case_id="adhoc")

    state = AgentState(case=manifest)

    planner: object | None = None
    if use_llm:
        from ..agent.config import load_config
        from ..agent.llm_factory import build_llm_client
        from ..agent.planner import LLMPlanner

        cfg = load_config()
        client = build_llm_client(cfg.llm)
        planner = LLMPlanner(client)  # ajusta si la firma requiere fallback
    else:
        planner = RuleBasedPlanner()

    ctx = SessionContext(state=state)
    return ctx, planner


def run_chat(
    case_path: Path | None = None,
    use_llm: bool = False,
    verbose: bool = False,
) -> int:
    ctx, planner = _build_context(case_path, use_llm)
    loop = AgentLoop(planner=planner)  # type: ignore[arg-type]

    SESSIONS_DIR.mkdir(exist_ok=True)
    sid = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    sdir = SESSIONS_DIR / sid
    sdir.mkdir(exist_ok=True)

    console.print(
        f"[dim]session {sid} · "
        f"planner={'llm' if use_llm else 'rules'} · "
        f"'exit' o Ctrl-C para salir[/dim]"
    )

    try:
        while True:
            try:
                user = console.input("[bold cyan]you ›[/bold cyan] ").strip()
            except (EOFError, KeyboardInterrupt):
                break
            if not user or user.lower() in {"exit", "quit", ":q"}:
                break

            ctx.state.note(f"user: {user}")

            try:
                for event in loop.stream(ctx):
                    render(event, console, verbose=verbose)
            except Exception as e:  # noqa: BLE001
                render(
                    Event(kind="error", payload={"message": str(e)}),
                    console,
                    verbose=verbose,
                )
    finally:
        out = sdir / "session.json"
        ctx.state.to_json(out)
        console.print(f"[dim]sesión guardada en {out}[/dim]")

    return 0
