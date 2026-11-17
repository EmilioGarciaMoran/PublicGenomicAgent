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
    """Construye un grafo local del ROI (vg autoindex: referencia + VCF) y sus índices de giraffe.

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

    # Un único grafo para todo: vg autoindex construye el grafo (referencia + VCF) y los
    # índices de giraffe; el .vg y el .xg se derivan del mismo GBZ. Así giraffe, la
    # proyección a BAM, vg pack y vg call usan la misma numeración de nodos.
    # (include_alt_paths y max_node_length ya no se usan: autoindex fija la construcción.)
    vcf_gz = vcf
    if not str(vcf).endswith(".gz"):
        vcf_gz = out_vg.with_suffix(".cohort.vcf.gz")
        runtime.run("bcftools", ["view", str(vcf), "-Oz", "-o", str(vcf_gz)])
    if not Path(str(vcf_gz) + ".tbi").exists() and not Path(str(vcf_gz) + ".csi").exists():
        runtime.run("bcftools", ["index", "-t", str(vcf_gz)])
    prefix = out_vg.with_suffix("")
    tmp = out_vg.parent / f".{out_vg.stem}_tmp"
    tmp.mkdir(exist_ok=True)
    runtime.run("vg", ["autoindex", "--workflow", "giraffe", "-r", str(ref), "-v", str(vcf_gz),
                       "-p", str(prefix), "-T", str(tmp)], stderr=subprocess.PIPE)
    gbz = Path(f"{prefix}.giraffe.gbz")
    out_xg = out_vg.with_suffix(".xg")
    with out_xg.open("wb") as f:
        runtime.run("vg", ["convert", "-x", str(gbz)], stdout=f, stderr=subprocess.PIPE)
    with out_vg.open("wb") as f:
        runtime.run("vg", ["convert", "-p", str(gbz)], stdout=f, stderr=subprocess.PIPE)

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
        outputs={"vg": str(out_vg), "xg": str(out_xg), "giraffe_gbz": str(gbz)},
        graph_vg=out_vg,
        graph_xg=out_xg,
        giraffe_gbz=gbz,
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
    """Alinea lecturas (FASTQ, simples o pareadas) contra el grafo local con vg giraffe.

    - output_bam: lecturas proyectadas (surject) a la referencia lineal, ordenadas e
      indexadas; es la entrada de call_from_graph(method="surject").
    - output_gam / output_pack: alineamiento en el grafo y su cobertura, solo para
      call_from_graph(method="graph_genotype").
    Requiere los índices de giraffe de build_local_graph (<prefijo>.giraffe.gbz).
    """
    _require_vg(runtime)
    g = Path(inp.graph_vg).expanduser().resolve()
    fq = Path(inp.reads_fastq).expanduser().resolve()
    fq2 = Path(inp.reads_fastq_2).expanduser().resolve() if inp.reads_fastq_2 else None
    if not g.exists():
        raise FileNotFoundError(f"grafo no existe: {g}")
    for f in (fq, fq2):
        if f is not None and not f.exists():
            raise FileNotFoundError(f"FASTQ no existe: {f}")
    if inp.output_bam is None and inp.output_gam is None:
        raise ValueError("indica output_bam (llamada por surject) y/o output_gam")
    if inp.output_pack is not None and inp.output_gam is None:
        raise ValueError("output_pack requiere output_gam")
    prefix = g.with_suffix("")
    idx = {"-Z": Path(f"{prefix}.giraffe.gbz"), "-m": Path(f"{prefix}.shortread.withzip.min"),
           "-z": Path(f"{prefix}.shortread.zipcodes"), "-d": Path(f"{prefix}.dist")}
    missing = [str(p) for p in idx.values() if not p.exists()]
    if missing:
        raise FileNotFoundError("faltan índices de giraffe (reconstruye con build_local_graph): "
                                + ", ".join(missing))
    base = [a for k, p in idx.items() for a in (k, str(p))] + ["-f", str(fq)]
    if fq2:
        base += ["-f", str(fq2)]

    bam: Path | None = None
    if inp.output_bam is not None:
        bam = Path(inp.output_bam).expanduser().resolve()
        bam.parent.mkdir(parents=True, exist_ok=True)
        raw = bam.with_suffix(".unsorted.bam")
        with raw.open("wb") as out_f:
            runtime.run("vg", ["giraffe", *base, "-o", "BAM", "--sample", inp.sample_name,
                               "--read-group", inp.sample_name], stdout=out_f, stderr=subprocess.PIPE)
        runtime.run("samtools", ["sort", "-o", str(bam), str(raw)])
        runtime.run("samtools", ["index", str(bam)])
        raw.unlink()

    gam: Path | None = None
    reads_total = reads_aligned = 0
    if inp.output_gam is not None:
        gam = Path(inp.output_gam).expanduser().resolve()
        gam.parent.mkdir(parents=True, exist_ok=True)
        with gam.open("wb") as out_f:
            runtime.run("vg", ["giraffe", *base, "-o", "gam"], stdout=out_f, stderr=subprocess.PIPE)
        try:
            for line in runtime.capture("vg", ["stats", "-a", str(gam)]).splitlines():
                key, _, val = line.partition(":")
                val = (val.strip().split() or ["0"])[0]
                if key.strip().lower() == "total alignments":
                    reads_total = int(val)
                elif key.strip().lower() == "total aligned":
                    reads_aligned = int(val)
        except Exception:  # noqa: BLE001
            pass
    elif bam is not None:
        try:
            st = runtime.capture("samtools", ["flagstat", str(bam)]).splitlines()
            reads_total = int(st[0].split()[0])
            reads_aligned = int(next(l for l in st if " mapped (" in l).split()[0])
        except Exception:  # noqa: BLE001
            pass

    pack_out: Path | None = None
    if inp.output_pack is not None:
        xg = g.with_suffix(".xg")
        if not xg.exists():
            runtime.run("vg", ["index", "-x", str(xg), "-L", str(g)])
        pack_out = Path(inp.output_pack).expanduser().resolve()
        pack_out.parent.mkdir(parents=True, exist_ok=True)
        runtime.run("vg", ["pack", "-x", str(xg), "-g", str(gam),
                           "-Q", str(inp.pack_min_quality), "-o", str(pack_out)])

    return AlignToGraphOutput(
        ok=True,
        tool="align_to_graph",
        message=f"{reads_aligned}/{reads_total} lecturas alineadas al grafo",
        outputs={k: str(v) for k, v in (("gam", gam), ("bam", bam), ("pack", pack_out)) if v},
        gam=gam,
        bam=bam,
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
    """Llama variantes tras alinear al grafo local.

    method="surject" (por defecto): bcftools mpileup + call sobre el BAM proyectado
    a la referencia lineal. Descubre variantes nuevas: la ganancia viene de alinear
    mejor, no de que el grafo contenga la respuesta (docs/benchmark_local_pangenome.md).

    method="graph_genotype": vg call sobre el pack. Solo genotipa sitios que ya están
    en el grafo (alelos de la cohorte); NO descubre variantes nuevas. Con un VCF de
    cohorte sin fase puede sesgar genotipos heterocigotos hacia 1/1.
    """
    vcf = Path(inp.output_vcf).expanduser().resolve()
    vcf.parent.mkdir(parents=True, exist_ok=True)
    vcf_plain = vcf.with_suffix(".plain.vcf")
    if inp.method == "surject":
        if inp.bam is None or inp.reference_fasta is None:
            raise ValueError("method='surject' requiere bam y reference_fasta")
        bam = Path(inp.bam).expanduser().resolve()
        ref = Path(inp.reference_fasta).expanduser().resolve()
        for f in (bam, ref):
            if not f.exists():
                raise FileNotFoundError(f"no existe: {f}")
        pile = vcf.with_suffix(".pileup.bcf")
        runtime.run("bcftools", ["mpileup", "-f", str(ref), "-a", "AD,DP", "-q", str(inp.min_mapq),
                                 "-d", "10000", "-Ou", "-o", str(pile), str(bam)])
        runtime.run("bcftools", ["call", "-mv", "-Ov", "-o", str(vcf_plain), str(pile)])
        pile.unlink()
    else:
        _require_vg(runtime)
        if not (inp.graph_vg and inp.graph_xg and inp.pack):
            raise ValueError("method='graph_genotype' requiere graph_vg, graph_xg y pack")
        g = Path(inp.graph_vg).expanduser().resolve()
        xg = Path(inp.graph_xg).expanduser().resolve()
        pack = Path(inp.pack).expanduser().resolve()
        for f in (g, xg, pack):
            if not f.exists():
                raise FileNotFoundError(f"no existe: {f}")
        raw = vcf.with_suffix(".vgcall.vcf")
        with raw.open("wb") as out_f:
            runtime.run("vg", ["call", str(xg), "-k", str(pack)],
                        stdout=out_f, stderr=subprocess.PIPE)
        # vg call informa por burbuja: variantes cercanas pueden salir en un solo
        # registro (REF de varias bases). Se separan en variantes simples.
        runtime.run("bcftools", ["norm", "-a", "-m", "-any", "-Ov", "-o", str(vcf_plain), str(raw)])
        raw.unlink()

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
        message=f"{variants_total} variantes ({inp.method})",
        outputs={"vcf": str(vcf_gz), "tbi": str(vcf_tbi)},
        vcf=vcf_gz,
        vcf_tbi=vcf_tbi,
        method=inp.method,
        variants_total=variants_total,
    )
