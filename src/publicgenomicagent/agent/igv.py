"""IGV (Integrative Genomics Viewer) integration.

Two mechanisms:

  1. Session XML: a self-contained XML file with reference
     and track paths that the user can open directly in IGV.
     Works for both IGV Desktop and igv.js.

  2. Control URLs: HTTP URLs to the local IGV Desktop
     port (default 60151) that trigger actions like goto,
     load, genome. Works only when IGV Desktop is running.

No external dependencies: string building only.
"""
from __future__ import annotations

import html as _html
import json
from pathlib import Path
from xml.sax.saxutils import escape as _xml_escape


DEFAULT_IGV_PORT = 60151


def _abs_uri(path: Path) -> str:
    """Devuelve un file:// URI absoluto para IGV."""
    return path.resolve().as_uri()


def build_control_url(
    *,
    locus: str | None = None,
    load_path: Path | None = None,
    genome: str | None = None,
    port: int = DEFAULT_IGV_PORT,
) -> str:
    """Construye una URL de control para IGV Desktop.

    IGV Desktop escucha en http://localhost:<port> y acepta
    parámetros GET. Ejemplos:
      /goto?locus=chr1:1000-1200
      /load?file=/path/to/file.bam
      /genome?id=hg38

    Solo una acción por URL.
    """
    base = f"http://localhost:{port}"
    if locus:
        return f"{base}/goto?locus={_html.escape(locus)}"
    if load_path is not None:
        return f"{base}/load?file={_html.escape(str(load_path.resolve()))}"
    if genome:
        return f"{base}/genome?id={_html.escape(genome)}"
    return base


def build_session_xml(
    *,
    reference_path: Path,
    tracks: list[tuple[str, Path]],
    locus: str | None = None,
) -> str:
    """Genera un session.xml de IGV con referencia y tracks.

    `tracks` es una lista de (nombre, path). Los paths se
    convierten a file:// URIs absolutos.

    Formato simplificado (compatible con IGV Desktop y igv.js):

      <Session genome="hg38" locus="chr1:1000-1200">
        <Resources>
          <Resource path="/abs/path/file.bam" name="..."/>
        </Resources>
      </Session>

    Para referencias, IGV espera un `genome` attribute, no un
    Resource. Cuando se especifica un reference_path, se añade
    como track adicional y se mantiene el genome="hg38" por
    defecto (asumiendo que la referencia coincide).
    """
    lines = []
    lines.append('<?xml version="1.0" encoding="UTF-8"?>')

    # Atributos del Session
    attrs = ['genome="hg38"']
    if locus:
        attrs.append(f'locus="{_xml_escape(locus)}"')
    lines.append(f'<Session {" ".join(attrs)}>')

    # Resources: tracks
    lines.append("  <Resources>")
    for name, path in tracks:
        uri = _abs_uri(Path(path))
        lines.append(
            f'    <Resource path="{_xml_escape(uri)}" '
            f'name="{_xml_escape(name)}" />'
        )
    lines.append("  </Resources>")

    # Panels: los tracks se muestran en orden
    lines.append("  <Panel name='Data Panel'>")
    for i, (name, _path) in enumerate(tracks):
        lines.append(f'    <Track id="{_xml_escape(name)}" />')
    lines.append("  </Panel>")

    lines.append("</Session>")
    return "\n".join(lines) + "\n"


def build_session_from_state(session_json: Path) -> dict:
    """Lee un session.json y construye todo lo necesario para IGV.

    Devuelve un dict con:
      - xml: string con el session.xml
      - locus: la primera ROI como string "chr:start-end"
      - goto_url: URL de control para saltar al locus
      - track_urls: lista de URLs de control para cada track
      - tracks: lista de (name, path)
    """
    data = json.loads(Path(session_json).read_text(encoding="utf-8"))
    case = data.get("case", {})
    bams = data.get("bams", {})
    vcfs = data.get("vcfs", {})

    # Locus: primera ROI
    locus = None
    rois = case.get("candidate_rois") or []
    if rois:
        roi = rois[0]
        locus = f"{roi.get('chrom')}:{roi.get('start')}-{roi.get('end')}"

    # Tracks: BAMs + VCFs
    tracks: list[tuple[str, Path]] = []
    for sample, path in bams.items():
        p = Path(path)
        if p.exists():
            tracks.append((f"{sample}.bam", p))
    for label, path in vcfs.items():
        p = Path(path)
        if p.exists():
            tracks.append((f"{label}.vcf.gz", p))

    # XML (necesitamos reference_path por compatibilidad, aunque no se use)
    refs = data.get("references", {})
    ref_path = next((Path(v) for v in refs.values() if Path(v).exists()), Path("."))
    xml = build_session_xml(
        reference_path=ref_path,
        tracks=tracks,
        locus=locus,
    )

    # URLs de control
    goto_url = build_control_url(locus=locus) if locus else build_control_url()
    track_urls = [
        {"name": name, "url": build_control_url(load_path=path)}
        for name, path in tracks
    ]

    return {
        "xml": xml,
        "locus": locus,
        "goto_url": goto_url,
        "track_urls": track_urls,
        "tracks": [(n, str(p)) for n, p in tracks],
    }

