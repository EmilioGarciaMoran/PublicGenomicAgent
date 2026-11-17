"""E3: qué muestra el visor a partir del session.json (sin herramientas externas)."""
from __future__ import annotations

from publicgenomicagent.agent.viewer import PAD, build_viewspec


def _session(tmp_path, ranking=True):
    hom = tmp_path / "hom.vcf"
    hom.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n"
                   "chr22\t500\t.\tC\tA\nchr22\t900\t.\tG\tT\n")
    bams = {s: str(tmp_path / f"{s}.bam") for s in ("P1", "FA", "MO")}
    calls = [{"tool_name": "fetch_roi", "ok": True, "input": {"bam_path": bams[s]},
              "output": {"output_bam": str(tmp_path / f"{s}.roi.bam")}} for s in ("P1", "MO", "FA")]
    out = {"mendelian_filter": {"auto_rec_hom_vcf": str(hom)}}
    if ranking:
        out["prioritize_variants"] = {"ranking": [
            {"rank": 1, "chrom": "chr22", "pos": 900, "ref": "G", "alt": "T"}]}
    return {"case": {"pedigree": [{"sample": s} for s in ("FA", "MO", "P1", "GF")],
                     "candidate_rois": [{"chrom": "chr22", "start": 100, "end": 2000}]},
            "bams": bams, "tool_calls": calls, "outputs": out,
            "vcfs": {"trio": "j.vcf.gz", "frequency": "f.vcf.gz"},
            "references": {"hg38": str(tmp_path / "ref.fa")}}


def test_viewspec_orders_samples_by_ped_and_centres_on_top_candidate(tmp_path):
    vs = build_viewspec(_session(tmp_path))
    assert vs.samples == ["FA", "MO", "P1"]                    # orden del PED, GF sin BAM fuera
    assert str(vs.family_vcf) == "f.vcf.gz"                    # el VCF con frecuencias manda
    assert (vs.chrom, vs.start, vs.end) == ("chr22", 900 - PAD, 900 + PAD)
    assert [h[1] for h in vs.highlights] == [900, 500]         # sin duplicar la rank 1
    assert vs.highlights[0][4] == "rank 1"


def test_viewspec_falls_back_to_roi_and_none_when_empty(tmp_path):
    s = _session(tmp_path, ranking=False)
    s["outputs"] = {}
    vs = build_viewspec(s)
    assert (vs.start, vs.end) == (100, 2000) and vs.highlights == []
    s["case"]["candidate_rois"] = []
    assert build_viewspec(s) is None


def test_redact_absolute_paths_keeps_basenames_and_scripts():
    from publicgenomicagent.agent.report import redact_absolute_paths
    doc = ('<td>/home/egarmo/pga_e2/results/e2.vcf.gz</td><a href="file:///home/u/x.bam">x</a>'
           ' 2/3 sites <script>var p="/home/keep/me.js"</script>')
    out = redact_absolute_paths(doc)
    assert "/home/egarmo" not in out and "<td>e2.vcf.gz</td>" in out
    assert "/home/u/" not in out and "2/3 sites" in out and "/home/keep/me.js" in out
