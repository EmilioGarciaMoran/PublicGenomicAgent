"""Las tools de PGA (build_local_graph, align_to_graph, call_from_graph) reproducen el benchmark:
descubren una variante nueva en un haplotipo divergente que el grafo NO contiene."""
from __future__ import annotations

import gzip
import random
import subprocess
import sys
from pathlib import Path

import pytest

from publicgenomicagent.env.micromamba import bin_path

pytestmark = pytest.mark.integration
sys.path.insert(0, str(Path(__file__).parents[2] / "scripts"))


def _calls(vcf):
    with gzip.open(vcf, "rt") as f:
        return {(int(c[1]), c[3], c[4]) for c in (l.split("\t") for l in f if not l.startswith("#"))}


def test_tools_discover_novel_variant_only_with_surject(tmp_path):
    if not (bin_path("pga-pangenome", "vg").exists() and bin_path("pga-hts", "wgsim").exists()):
        pytest.skip("vg o wgsim no disponibles")
    import benchmark_local_pangenome as bm
    from publicgenomicagent.cli.main import _runtime
    from publicgenomicagent.tools.base import AlignToGraphInput, BuildLocalGraphInput, CallFromGraphInput
    from publicgenomicagent.tools.local_graph import align_to_graph, build_local_graph, call_from_graph

    ref, cohort, alts, _ = bm.design(random.Random(1), spacing=12)
    ref_fa, samp_fa, vcf = tmp_path / "ref.fa", tmp_path / "s.fa", tmp_path / "cohort.vcf"
    bm.write_fa(ref_fa, bm.CHROM, ref)
    bm.write_fa(samp_fa, bm.CHROM, bm.sample_seq(ref, cohort, alts))
    bm.write_vcf(vcf, cohort)                                   # sin las variantes objetivo
    subprocess.run([str(bin_path("pga-hts", "samtools")), "faidx", str(ref_fa)], check=True)
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    subprocess.run([str(bin_path("pga-hts", "wgsim")), "-S", "1", "-N", "3000", "-1", "150", "-2", "150",
                    "-d", "400", "-s", "40", "-e", "0.002", "-r", "0", "-R", "0",
                    str(samp_fa), str(r1), str(r2)], check=True, capture_output=True)
    rt = _runtime()
    g = build_local_graph(rt, BuildLocalGraphInput(reference_fasta=ref_fa, cohort_vcf=vcf,
                                                    output_vg=tmp_path / "g" / "local.vg"))
    assert g.giraffe_gbz and g.giraffe_gbz.exists()
    a = align_to_graph(rt, AlignToGraphInput(graph_vg=g.graph_vg, reads_fastq=r1, reads_fastq_2=r2,
                                             output_bam=tmp_path / "s.bam", output_gam=tmp_path / "s.gam",
                                             output_pack=tmp_path / "s.pack", sample_name="P1"))
    assert a.bam.exists() and a.reads_aligned > 0
    sj = call_from_graph(rt, CallFromGraphInput(output_vcf=tmp_path / "sj.vcf.gz", bam=a.bam,
                                                reference_fasta=ref_fa, min_mapq=20, sample_name="P1"))
    gg = call_from_graph(rt, CallFromGraphInput(output_vcf=tmp_path / "gg.vcf.gz", method="graph_genotype",
                                                graph_vg=g.graph_vg, graph_xg=g.graph_xg, pack=a.pack))
    t, r, alt = alts["divergent"]
    target = (t + 1, r, alt)
    assert target in _calls(sj.vcf)            # descubierta desde las lecturas
    assert target not in _calls(gg.vcf)        # vg call no puede: no está en el grafo
    assert sj.method == "surject" and gg.method == "graph_genotype"
