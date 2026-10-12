"""REPL conversacional sobre AgentLoop.stream()."""
from __future__ import annotations

import json
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

_REF_GLOBS = ("*.fa", "*.fasta", "*.fa.gz", "*.fasta.gz")


def _discover_bams(case_dir: Path) -> dict[str, Path]:
    """Convención: <case_dir>/bams/<sample>.bam -> {sample: path}."""
    bams_dir = case_dir / "bams"
    if not bams_dir.is_dir():
        return {}
    return {bam.stem: bam for bam in sorted(bams_dir.glob("*.bam"))}


def _discover_reference(case_dir: Path) -> dict[str, Path]:
    """Convención: <case_dir>/ref/<algo>.fa -> {"ref": path}."""
    ref_dir = case_dir / "ref"
    if not ref_dir.is_dir():
        return {}
    for pattern in _REF_GLOBS:
        for fa in sorted(ref_dir.glob(pattern)):
            return {"ref": fa}
    return {}


def _build_context(
    case_path: Path | None,
    use_llm: bool,
) -> tuple[SessionContext, object | None]:
    """Construye SessionContext + planner.

    Sin --case, se usa un CaseManifest mínimo ("adhoc") sin BAMs ni ROIs.
    Con --case, se cargan BAMs y referencia por convención de directorio
    (<case_dir>/bams/*.bam y <case_dir>/ref/*.fa*).
    Sin --llm, el planner es RuleBasedPlanner (determinista, sin red).
    """
    # Runtime bioinformático: resuelve binarios (samtools, bcftools, ...)
    # en sus entornos aislados. Reutilizamos el patrón de cli/main.py.
    from ..cli.main import _runtime as _build_runtime

    if case_path is not None:
        if not case_path.exists():
            raise FileNotFoundError(f"no existe el manifest: {case_path}")
        raw = json.loads(case_path.read_text(encoding="utf-8"))
        manifest = CaseManifest.model_validate(raw)

        case_dir = case_path.parent
        bams = _discover_bams(case_dir)
        references = _discover_reference(case_dir)

        state = AgentState(case=manifest, bams=bams, references=references)

        if bams:
            console.print(
                f"[dim]· {len(bams)} BAM(s) cargados: "
                f"{', '.join(sorted(bams))}[/dim]"
            )
        if references:
            console.print(
                f"[dim]· referencia: {next(iter(references.values()))}[/dim]"
            )
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
        planner = LLMPlanner(client)
    else:
        planner = RuleBasedPlanner()

    ctx = SessionContext(state=state, runtime=_build_runtime())
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
