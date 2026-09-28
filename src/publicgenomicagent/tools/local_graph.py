from __future__ import annotations

import re
import subprocess
from pathlib import Path

from ..env.runtime import ToolRuntime
from .base import (
    AlignToGraphInput,
    AlignToGraphOutput,
    BuildLocalGraphInput,
    BuildLocalGraphOutput,
    CallFromGraphInput,
    CallFromGraphOutput,
)


# ---------------------------------------------------------------------
# Utilidades comunes
# ---------------------------------------------------------------------

def _require_vg(runtime: ToolRuntime) -> None:
    """Verifica que vg está accesible vía ToolRuntime."""
    runtime.spec("vg")


def _count_gam_reads(gam: Path) -> int:
    """Cuenta lecturas en un GAM sin cargarlo todo en memoria.

    El GAM es binario; usamos `vg view -c` o contamos marcadores.
    Simplificación: usamos el tamaño como proxy y no contamos exacto.
    """
    # Implementación real requiere vg; se hace en align_to_graph
    return 0


def _count_vcf_variants(vcf: Path) -> int:
    """Cuenta variantes no-header en un VCF."""
    if not vcf.exists():
        return 0
    count = 0
    import gzip
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt") as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                count += 1
    return count


# ---------------------------------------------------------------------
# 1. build_local_graph
# ---------------------------------------------------------------------

def build_local_graph(
    runtime: ToolRuntime,
    inp: BuildLocalGraphInput,
) -> BuildLocalGraphOutput:
    """Construye un grafo local del ROI con vg construct.

    - reference_fasta: FASTA del ROI (con contig coherente con el VCF).
    - cohort_vcf: VCF con variantes poblacionales (AF, SVs, etc.).
    - output_vg: fichero .vg resultante.
    - include_alt_paths: si True, usa -a (preserva alt paths, necesario
      para vg call).
    """
    _require_vg(runtime)

    ref = Path(inp.reference_fasta).expanduser().resolve()
    vcf = Path(inp.cohort_vcf).expanduser().resolve()
    out_vg = Path(inp.output_vg).expanduser().resolve()

    if not ref.exists():
        raise FileNotFoundError(f"FASTA de referencia no existe: {ref}")
    if not vcf.exists():
        raise FileNotFoundError(f"cohort VCF no existe: {vcf}")

    out_vg.parent.mkdir(parents=True, exist_ok=True)

    # Comando base
    args: list[str] = [
        "construct",
        "-r", str(ref),
        "-v", str(vcf),
        "-m", str(inp.max_node_length),
    ]
    if inp.include_alt_paths:
        args.append("-a")

    # vg construct escribe a stdout; redirigimos
    with out_vg.open("wb") as out_f:
        result = runtime.run("vg", args, stdout=out_f, stderr=subprocess.PIPE)

    # vg construct devuelve un .vg; ahora indexamos
    out_xg = out_vg.with_suffix(".xg")
    index_args = ["index", "-x", str(out_xg), "-L", str(out_vg)]
    runtime.run("vg", index_args)

    # Contar nodos y edges con vg stats -z (devuelve dos líneas: nodes/edges)
    nodes, edges = 0, 0
    try:
        stats_out = runtime.capture("vg", ["stats", "-z", str(out_vg)])
        for line in stats_out.splitlines():
            parts = line.strip().split()
            if len(parts) >= 2:
                if parts[0] == "nodes":
                    nodes = int(parts[1])
                elif parts[0] == "edges":
                    edges = int(parts[1])
    except Exception:
        pass

    return BuildLocalGraphOutput(
        ok=True,
        tool="build_local_graph",
        message=f"Grafo {out_vg.name}: {nodes} nodos, {edges} edges",
        outputs={"vg": str(out_vg), "xg": str(out_xg)},
        graph_vg=out_vg,
        graph_xg=out_xg,
        nodes=nodes,
        edges=edges,
    )


# ---------------------------------------------------------------------
# 2. align_to_graph
# ---------------------------------------------------------------------

def align_to_graph(
    runtime: ToolRuntime,
    inp: AlignToGraphInput,
) -> AlignToGraphOutput:
    """Alinea lecturas (FASTQ) contra un grafo con vg giraffe.

    - graph_vg: grafo .vg.
    - reads_fastq: lecturas sin alinear.
    - output_gam: alineamiento en formato GAM.
    - output_pack: opcional, cobertura empaquetada para vg call.
    """
    _require_vg(runtime)

    g = Path(inp.graph_vg).expanduser().resolve()
    fq = Path(inp.reads_fastq).expanduser().resolve()
    gam = Path(inp.output_gam).expanduser().resolve()

    if not g.exists():
        raise FileNotFoundError(f"grafo no existe: {g}")
    if not fq.exists():
        raise FileNotFoundError(f"FASTQ no existe: {fq}")

    gam.parent.mkdir(parents=True, exist_ok=True)

    # giraffe puede crear automáticamente los índices secundarios
    # (gbz, dist, minimizers) si no están.
    xg = g.with_suffix(".xg")
    if not xg.exists():
        runtime.run("vg", ["index", "-x", str(xg), "-L", str(g)])

    # vg giraffe: lee FASTQ (-f), escribe GAM a stdout
    with gam.open("wb") as out_f:
        runtime.run(
            "vg",
            ["giraffe", "-x", str(xg), "-f", str(fq), "-o", "gam"],
            stdout=out_f,
            stderr=subprocess.PIPE,
        )

    # Contar lecturas alineadas con vg stats -a (o con vg view -a)
    reads_total = 0
    reads_aligned = 0
    try:
        stats_out = runtime.capture("vg", ["stats", "-a", str(gam)])
        for line in stats_out.splitlines():
            if ":" not in line:
                continue
            key, _, val = line.partition(":")
            key = key.strip().lower()
            val = val.strip().split()[0] if val.strip() else "0"
            if key == "total alignments":
                reads_total = int(val)
            elif key == "total aligned":
                reads_aligned = int(val)
    except Exception:
        pass

    # Pack (opcional)
    pack_out: Path | None = None
    if inp.output_pack is not None:
        pack_out = Path(inp.output_pack).expanduser().resolve()
        pack_out.parent.mkdir(parents=True, exist_ok=True)
        runtime.run(
            "vg",
            ["pack", "-x", str(xg), "-g", str(gam),
             "-Q", str(inp.pack_min_quality), "-o", str(pack_out)],
        )

    return AlignToGraphOutput(
        ok=True,
        tool="align_to_graph",
        message=f"{reads_aligned}/{reads_total} lecturas alineadas al grafo",
        outputs={
            "gam": str(gam),
            **({"pack": str(pack_out)} if pack_out else {}),
        },
        gam=gam,
        pack=pack_out,
        reads_aligned=reads_aligned,
        reads_total=reads_total,
    )


# ---------------------------------------------------------------------
# 3. call_from_graph
# ---------------------------------------------------------------------

def call_from_graph(
    runtime: ToolRuntime,
    inp: CallFromGraphInput,
) -> CallFromGraphOutput:
    """Llama variantes desde un pack sobre un grafo con vg call.

    - graph_vg: grafo .vg.
    - graph_xg: índice .xg.
    - pack: cobertura empaquetada (de align_to_graph).
    - output_vcf: VCF resultante.
    """
    _require_vg(runtime)

    g = Path(inp.graph_vg).expanduser().resolve()
    xg = Path(inp.graph_xg).expanduser().resolve()
    pack = Path(inp.pack).expanduser().resolve()
    vcf = Path(inp.output_vcf).expanduser().resolve()

    if not g.exists():
        raise FileNotFoundError(f"grafo no existe: {g}")
    if not xg.exists():
        raise FileNotFoundError(f"índice .xg no existe: {xg}")
    if not pack.exists():
        raise FileNotFoundError(f"pack no existe: {pack}")

    vcf.parent.mkdir(parents=True, exist_ok=True)

    # vg call emite VCF plano por stdout
    vcf_plain = vcf.with_suffix(".plain.vcf")
    with vcf_plain.open("wb") as out_f:
        runtime.run(
            "vg",
            ["call", str(xg), "-k", str(pack)],
            stdout=out_f,
            stderr=subprocess.PIPE,
        )

    # Comprimir a .vcf.gz con bcftools (BGZF) para permitir tabix
    # El nombre final esperado es <...>.vcf.gz; aceptamos tanto .vcf
    # como .vcf.gz en el input y normalizamos.
    vcf_gz = vcf if str(vcf).endswith(".gz") else vcf.with_suffix(vcf.suffix + ".gz")
    runtime.run(
        "bcftools",
        ["view", str(vcf_plain), "-Oz", "-o", str(vcf_gz)],
    )
    vcf_plain.unlink()

    # Renombrar el sample si se ha proporcionado un nombre
    # (vg call produce "SAMPLE" por defecto)
    if inp.sample_name:
        # bcftools reheader -s <(echo NAME) es la forma limpia,
        # pero desde Python es más simple escribir un fichero temporal
        # con el nombre y usar reheader -s.
        samples_file = vcf_gz.parent / "_samples.txt"
        samples_file.write_text(inp.sample_name + "\n")
        vcf_gz_renamed = vcf_gz.with_name(vcf_gz.stem + "_renamed.vcf.gz")
        runtime.run(
            "bcftools",
            ["reheader", "-s", str(samples_file),
             str(vcf_gz), "-o", str(vcf_gz_renamed)],
        )
        vcf_gz.unlink()
        vcf_gz_renamed.rename(vcf_gz)
        samples_file.unlink()

    # Indexar
    vcf_tbi = Path(str(vcf_gz) + ".tbi")
    runtime.run("bcftools", ["index", "-t", str(vcf_gz)])

    variants_total = _count_vcf_variants(vcf_gz)

    return CallFromGraphOutput(
        ok=True,
        tool="call_from_graph",
        message=f"{variants_total} variantes llamadas desde el grafo",
        outputs={"vcf": str(vcf_gz), "tbi": str(vcf_tbi)},
        vcf=vcf_gz,
        vcf_tbi=vcf_tbi,
        variants_total=variants_total,
    )
