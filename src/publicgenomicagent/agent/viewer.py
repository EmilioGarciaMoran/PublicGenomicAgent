"""Visor embebido del informe: ViewSpec neutra + adaptador JBrowse autocontenido.

`build_viewspec(session)` decide QUÉ se muestra a partir del session.json
(locus, referencia, sub-BAM por muestra en orden del PED, VCF familiar y
candidatas). Es independiente del visor.

`render_jbrowse(viewspec, workdir)` recorta los datos a una ventana
alrededor del locus y devuelve el HTML del visor JBrowse con todo dentro
(componente MIT/Apache compilado en vendor/, datos en base64): abre sin
red con file://. La referencia se incrusta como el cromosoma completo con
su longitud real, secuencia real en la ventana y N en el resto, comprimido
con bgzip: las coordenadas son las reales y pesa unos cientos de KB.

Aviso de privacidad: si se incluyen lecturas (include_reads=True), el
informe contiene datos genómicos individuales; compártase solo con casos
sintéticos o por un canal controlado.
"""
from __future__ import annotations

import base64
import gzip
import hashlib
import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from ..env.micromamba import bin_path

VENDOR = Path(__file__).parent / "vendor"
BUNDLE = VENDOR / "jbrowse-lgv.min.js"
PAD = 150            # pb a cada lado de la variante para el locus inicial
WINDOW = 5_000       # pb a cada lado del locus que se incrustan


@dataclass
class ViewSpec:
    chrom: str
    start: int
    end: int
    reference: Path | None
    samples: list[str] = field(default_factory=list)          # orden del PED
    sub_bams: dict[str, Path] = field(default_factory=dict)   # sample -> sub-BAM de la ROI
    family_vcf: Path | None = None
    highlights: list[tuple[str, int, str, str, str]] = field(default_factory=list)  # chrom,pos,ref,alt,motivo

    @property
    def locus(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


def _vcf_records(path: Path, limit: int = 200) -> list[tuple[str, int, str, str]]:
    out = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            out.append((c[0], int(c[1]), c[3], c[4]))
            if len(out) >= limit:
                break
    return out


def build_viewspec(session: dict) -> ViewSpec | None:
    """Vista a partir del session.json. None si no hay nada que mostrar."""
    case = session.get("case", {})
    bams = {k: str(v) for k, v in (session.get("bams") or {}).items()}
    by_path = {v: k for k, v in bams.items()}
    sub = {}
    for c in session.get("tool_calls", []):
        if c.get("tool_name") == "fetch_roi" and c.get("ok") and c.get("output"):
            s = by_path.get(c["input"].get("bam_path"))
            if s and c["output"].get("output_bam"):
                sub[s] = Path(c["output"]["output_bam"])
    order = [i["sample"] for i in case.get("pedigree", []) if i.get("sample") in sub] or list(sub)
    vcfs = session.get("vcfs") or {}
    fam_vcf = next((Path(vcfs[k]) for k in ("frequency", "annotated", "trio") if vcfs.get(k)), None)

    highlights: list[tuple[str, int, str, str, str]] = []
    outputs = session.get("outputs") or {}
    for rk in (outputs.get("prioritize_variants") or {}).get("ranking", [])[:10]:
        highlights.append((rk["chrom"], int(rk["pos"]), rk["ref"], rk["alt"], f"rank {rk['rank']}"))
    mf = outputs.get("mendelian_filter") or {}
    for key, why in (("auto_rec_hom_vcf", "recessive homozygous"), ("comp_het_vcf", "compound heterozygous")):
        p = mf.get(key)
        if p and Path(p).exists():
            for chrom, pos, ref, alt in _vcf_records(Path(p)):
                if not any(h[1] == pos and h[3] == alt for h in highlights):
                    highlights.append((chrom, pos, ref, alt, why))

    rois = case.get("candidate_rois") or []
    if highlights:
        chrom, pos = highlights[0][0], highlights[0][1]
        start, end = max(1, pos - PAD), pos + PAD
    elif rois:
        chrom, start, end = rois[0]["chrom"], int(rois[0]["start"]), int(rois[0]["end"])
    else:
        return None
    ref = (session.get("references") or {}).get("hg38") or (session.get("references") or {}).get("ref")
    return ViewSpec(chrom=chrom, start=start, end=end, reference=Path(ref) if ref else None,
                    samples=order, sub_bams=sub, family_vcf=fam_vcf, highlights=highlights)


def _run(cmd: list, **kw) -> subprocess.CompletedProcess:
    return subprocess.run([str(c) for c in cmd], check=True, capture_output=True, **kw)


def _b64(p: Path) -> str:
    return base64.b64encode(Path(p).read_bytes()).decode()


def _window_reference(vs: ViewSpec, ws: int, we: int, out: Path) -> Path:
    """Cromosoma de longitud real: secuencia real en [ws, we], N en el resto; bgzip + fai + gzi."""
    samtools, bgzip = bin_path("pga-hts", "samtools"), bin_path("pga-hts", "bgzip")
    fai = Path(str(vs.reference) + ".fai")
    length = next(int(l.split("\t")[1]) for l in fai.read_text().splitlines()
                  if l.split("\t")[0] == vs.chrom)
    seq = "".join(_run([samtools, "faidx", vs.reference, f"{vs.chrom}:{ws}-{we}"],
                       text=True).stdout.splitlines()[1:])
    plain = out.with_suffix("")                       # .fa
    with plain.open("w") as f:
        f.write(f">{vs.chrom}\n")
        col = 0

        def put(s: str):
            nonlocal col
            i = 0
            while i < len(s):
                take = min(len(s) - i, 60 - col)
                f.write(s[i:i + take]); col += take; i += take
                if col == 60:
                    f.write("\n"); col = 0
        for n in (ws - 1, None, length - we):
            if n is None:
                put(seq)
                continue
            while n:
                chunk = min(n, 600_000)
                put("N" * chunk); n -= chunk
        if col:
            f.write("\n")
    _run([bgzip, "-f", "-i", plain])                    # -> .fa.gz + .fa.gz.gzi
    _run([samtools, "faidx", out])                      # -> .fa.gz.fai
    return out


def render_jbrowse(vs: ViewSpec, workdir: Path, include_reads: bool = True) -> tuple[str, dict]:
    """HTML (div + datos + script) del visor y metadatos (tamaños, hash del componente)."""
    if not BUNDLE.exists():
        raise FileNotFoundError(f"componente JBrowse no encontrado: {BUNDLE}")
    if vs.reference is None or not Path(vs.reference).exists():
        raise FileNotFoundError("referencia de la sesión no disponible para el visor")
    workdir.mkdir(parents=True, exist_ok=True)
    samtools, bcftools = bin_path("pga-hts", "samtools"), bin_path("pga-hts", "bcftools")
    ws, we = max(1, vs.start - WINDOW), vs.end + WINDOW
    region = f"{vs.chrom}:{ws}-{we}"
    files: dict[str, Path] = {}

    ref = _window_reference(vs, ws, we, workdir / "ref_window.fa.gz")
    files.update({"ref.fa.gz": ref, "ref.fa.gz.fai": Path(str(ref) + ".fai"),
                  "ref.fa.gz.gzi": Path(str(ref) + ".gzi")})

    def loc(name: str) -> dict:
        return {"locationType": "BlobLocation", "blobId": name, "name": name}

    tracks, show = [], []
    if vs.family_vcf and Path(vs.family_vcf).exists():
        v = workdir / "family_window.vcf.gz"
        _run([bcftools, "view", "-r", region, "-Oz", "-o", v, vs.family_vcf])
        _run([bcftools, "index", "-f", "-t", v])
        files.update({"family.vcf.gz": v, "family.vcf.gz.tbi": Path(str(v) + ".tbi")})
        tracks.append({"type": "VariantTrack", "trackId": "family_vcf", "name": "Family genotypes",
                       "assemblyNames": ["GRCh38"],
                       "adapter": {"type": "VcfTabixAdapter", "vcfGzLocation": loc("family.vcf.gz"),
                                   "index": {"location": loc("family.vcf.gz.tbi")}},
                       "displays": [{"type": "MultiLinearVariantDisplay",
                                     "displayId": "family_vcf-multi"}]})
        show.append("family_vcf")
    if include_reads:
        for s in vs.samples:
            b = workdir / f"{s}.window.bam"
            _run([samtools, "view", "-b", "-o", b, vs.sub_bams[s], region])
            _run([samtools, "index", b])
            files.update({f"{s}.bam": b, f"{s}.bam.bai": Path(str(b) + ".bai")})
            tracks.append({"type": "AlignmentsTrack", "trackId": f"bam_{s}", "name": f"{s} (reads)",
                           "assemblyNames": ["GRCh38"],
                           "adapter": {"type": "BamAdapter", "bamLocation": loc(f"{s}.bam"),
                                       "index": {"location": loc(f"{s}.bam.bai")}}})
            show.append(f"bam_{s}")

    data = {
        "assembly": {"name": "GRCh38", "sequence": {
            "type": "ReferenceSequenceTrack", "trackId": "refseq",
            "adapter": {"type": "BgzipFastaAdapter", "fastaLocation": loc("ref.fa.gz"),
                        "faiLocation": loc("ref.fa.gz.fai"), "gziLocation": loc("ref.fa.gz.gzi")}}},
        "tracks": tracks, "show": show, "locus": vs.locus,
        "files": {k: _b64(p) for k, p in files.items()},
    }
    bundle = BUNDLE.read_text(encoding="utf-8")
    payload = json.dumps(data).replace("</", "<\\/")
    html = ('<div id="pga-jbrowse" style="min-height:420px"></div>'
            f'<script id="pga-jbrowse-data" type="application/json">{payload}</script>'
            f"<script>{bundle}</script>")
    meta = {
        "locus": vs.locus, "window": region, "include_reads": include_reads,
        "embedded_bytes": sum(Path(p).stat().st_size for p in files.values()),
        "bundle_sha256": hashlib.sha256(BUNDLE.read_bytes()).hexdigest(),
        "tracks": [t["name"] for t in tracks],
    }
    return html, meta
