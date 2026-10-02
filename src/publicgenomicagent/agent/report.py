"""HTML report generator for a PublicGenomicAgent session.

Reads a session.json (as written by AgentState.to_json) and
produces a self-contained HTML file with:

  - Header: case_id, source, proband, generation timestamp.
  - Case summary: pedigree, ROIs, HPO terms.
  - Pipeline: table of tool calls with status and duration.
  - Results: VCFs and their logical keys.
  - Candidate variants: counts from mendelian_filter.
  - Privacy footer.

No external dependencies. No Jinja. Inline CSS only.
"""
from __future__ import annotations

import html as _html
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


CSS = """
  :root { --fg: #1a202c; --bg: #f7fafc; --accent: #2b6cb0; --ok: #38a169; --err: #e53e3e; --muted: #718096; --border: #e2e8f0; }
  * { box-sizing: border-box; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 0; padding: 2rem; background: var(--bg); color: var(--fg); line-height: 1.5; }
  .container { max-width: 960px; margin: 0 auto; background: white; border-radius: 8px; padding: 2rem 2.5rem; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }
  h1 { font-size: 1.75rem; margin: 0 0 0.25rem 0; }
  h2 { font-size: 1.15rem; margin: 2rem 0 0.75rem 0; padding-bottom: 0.4rem; border-bottom: 1px solid var(--border); }
  .meta { color: var(--muted); font-size: 0.875rem; margin-bottom: 1.5rem; }
  table { width: 100%; border-collapse: collapse; margin: 0.5rem 0; font-size: 0.9rem; }
  th { text-align: left; padding: 0.5rem 0.75rem; background: var(--bg); border-bottom: 2px solid var(--border); font-weight: 600; }
  td { padding: 0.5rem 0.75rem; border-bottom: 1px solid var(--border); }
  code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; background: var(--bg); padding: 0.1em 0.35em; border-radius: 3px; font-size: 0.875em; }
  .ok { color: var(--ok); font-weight: 600; }
  .err { color: var(--err); font-weight: 600; }
  .badge { display: inline-block; padding: 0.15rem 0.5rem; border-radius: 4px; font-size: 0.75rem; font-weight: 600; }
  .badge.ok { background: #c6f6d5; color: #22543d; }
  .badge.err { background: #fed7d7; color: #742a2a; }
  .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 0.75rem; margin: 0.75rem 0; }
  .stat { padding: 0.75rem 1rem; background: var(--bg); border-radius: 6px; }
  .stat .label { font-size: 0.75rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.05em; }
  .stat .value { font-size: 1.5rem; font-weight: 700; margin-top: 0.25rem; }
  .privacy { margin-top: 2.5rem; padding: 1rem; background: #f0fff4; border-left: 3px solid var(--ok); font-size: 0.875rem; color: #22543d; }
  .empty { color: var(--muted); font-style: italic; padding: 0.5rem 0; }
  a { color: var(--accent); text-decoration: none; }
  a:hover { text-decoration: underline; }
"""


def _e(value: Any) -> str:
    """Escape HTML, tolerating non-string types."""
    if value is None:
        return "—"
    return _html.escape(str(value))


def _row(cells: list[str]) -> str:
    tds = "".join(f"<td>{c}</td>" for c in cells)
    return f"<tr>{tds}</tr>"


def _section(title: str, body: str) -> str:
    return f"<h2>{_e(title)}</h2>\n{body}\n"


def _render_pedigree(case: dict) -> str:
    pedigree = case.get("pedigree") or []
    if not pedigree:
        return '<p class="empty">No pedigree provided.</p>'
    rows = []
    for ind in pedigree:
        rel = []
        if ind.get("father"):
            rel.append(f"father={_e(ind['father'])}")
        if ind.get("mother"):
            rel.append(f"mother={_e(ind['mother'])}")
        status = "affected" if ind.get("affected") else "unaffected"
        rows.append(_row([
            f"<code>{_e(ind.get('sample'))}</code>",
            _e(ind.get("sex", "U")),
            status,
            ", ".join(rel) if rel else "—",
        ]))
    return (
        "<table><thead><tr><th>Sample</th><th>Sex</th>"
        "<th>Status</th><th>Relations</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def _render_rois(case: dict) -> str:
    rois = case.get("candidate_rois") or []
    if not rois:
        return '<p class="empty">No ROIs provided.</p>'
    rows = []
    for roi in rois:
        label = roi.get("label") or "—"
        loc = f"{roi.get('chrom')}:{roi.get('start')}-{roi.get('end')}"
        rows.append(_row([f"<code>{_e(loc)}</code>", _e(label)]))
    return (
        "<table><thead><tr><th>Region</th><th>Label</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def _render_tool_calls(tool_calls: list[dict]) -> str:
    if not tool_calls:
        return '<p class="empty">No tool calls recorded.</p>'
    rows = []
    for i, c in enumerate(tool_calls):
        status = (
            '<span class="badge ok">ok</span>'
            if c.get("ok")
            else f'<span class="badge err">error</span>'
        )
        err = c.get("error") or "—"
        duration = c.get("duration_ms", 0)
        rows.append(_row([
            str(i),
            f"<code>{_e(c.get('tool_name'))}</code>",
            status,
            f"{duration} ms",
            f"<span style='color:var(--muted);font-size:0.8em'>{_e(err)}</span>" if err != "—" else "—",
        ]))
    return (
        "<table><thead><tr><th>#</th><th>Tool</th><th>Status</th>"
        "<th>Duration</th><th>Error</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def _render_vcfs(vcfs: dict) -> str:
    if not vcfs:
        return '<p class="empty">No VCFs produced.</p>'
    rows = []
    for label, path in vcfs.items():
        p = Path(path)
        link = f"file://{p.resolve()}" if p.exists() else ""
        link_html = (
            f'<a href="{_e(link)}">open</a>' if link
            else '<span style="color:var(--muted)">missing</span>'
        )
        rows.append(_row([
            f"<code>{_e(label)}</code>",
            f"<code style='font-size:0.85em'>{_e(p.name)}</code>",
            link_html,
        ]))
    return (
        "<table><thead><tr><th>Key</th><th>File</th><th>Link</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def _render_mendelian_stats(outputs: dict) -> str:
    mf = outputs.get("mendelian_filter")
    if not mf:
        return '<p class="empty">mendelian_filter not executed.</p>'
    counts = mf.get("counts") or {}
    if not counts:
        return '<p class="empty">No counts available.</p>'
    stats = []
    for key in ["de_novo", "auto_rec_hom", "auto_dom", "x_linked_rec"]:
        v = counts.get(key, 0)
        stats.append(
            f'<div class="stat"><div class="label">{_e(key)}</div>'
            f'<div class="value">{_e(v)}</div></div>'
        )
    return f'<div class="grid">{"".join(stats)}</div>'


def render_session_report(session_json: Path | str, output_html: Path | str) -> Path:
    """Render a session.json into a self-contained HTML report.

    Returns the path to the written HTML file.
    """
    session_path = Path(session_json)
    output_path = Path(output_html)

    if not session_path.exists():
        raise FileNotFoundError(f"session.json not found: {session_path}")

    data = json.loads(session_path.read_text(encoding="utf-8"))
    case = data.get("case", {})
    tool_calls = data.get("tool_calls", [])
    vcfs = data.get("vcfs", {})
    outputs = data.get("outputs", {})

    case_id = case.get("case_id", "unknown")
    source = case.get("source", "—")
    proband = case.get("proband") or "—"
    n_ok = sum(1 for c in tool_calls if c.get("ok"))
    n_err = len(tool_calls) - n_ok

    gen_ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    # Build the HTML body
    body_parts = []
    body_parts.append("<h1>PublicGenomicAgent — Session Report</h1>")
    body_parts.append(
        f'<div class="meta">Case <code>{_e(case_id)}</code> · '
        f'Source {_e(source)} · Proband {_e(proband)} · '
        f'Generated {_e(gen_ts)}</div>'
    )
    body_parts.append(
        '<div class="grid">'
        f'<div class="stat"><div class="label">Tool calls</div><div class="value">{len(tool_calls)}</div></div>'
        f'<div class="stat"><div class="label">Success</div><div class="value" style="color:var(--ok)">{n_ok}</div></div>'
        f'<div class="stat"><div class="label">Errors</div><div class="value" style="color:var(--err)">{n_err}</div></div>'
        f'<div class="stat"><div class="label">VCFs</div><div class="value">{len(vcfs)}</div></div>'
        '</div>'
    )

    body_parts.append(_section("Pedigree", _render_pedigree(case)))
    body_parts.append(_section("Regions of interest", _render_rois(case)))
    body_parts.append(_section("Pipeline executed", _render_tool_calls(tool_calls)))
    body_parts.append(_section("Candidate variants", _render_mendelian_stats(outputs)))
    body_parts.append(_section("VCF files produced", _render_vcfs(vcfs)))

    body_parts.append(
        '<div class="privacy">'
        'Generated locally by PublicGenomicAgent. No network calls, '
        'no data left this machine. See <code>docs/privacy.md</code>.'
        '</div>'
    )

    document = (
        "<!DOCTYPE html>\n"
        "<html lang=\"en\"><head>"
        '<meta charset="utf-8">'
        f'<title>Report — {_e(case_id)}</title>'
        f"<style>{CSS}</style>"
        '</head><body>'
        '<div class="container">'
        + "".join(body_parts)
        + '</div></body></html>\n'
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(document, encoding="utf-8")
    return output_path

