#!/usr/bin/env python3
"""Genera fixture de cohorte: VCF bgzip + tabix + .fam + ground truth.

Coordenadas LOCALES al ROI (1-based, dentro del FASTA). El contig se
renombra a `chr2_roi` para evitar confusión con chr2 genómico.
"""
from __future__ import annotations
import json, random, subprocess, sys
from pathlib import Path

# NOTA: este script requiere el entorno pga-core.
# Ejecutar con:  ~/.pga/envs/pga-core/bin/python make_cohort_fixture.py <outdir>
# NO usar el python3 del sistema.
try:
    import publicgenomicagent  # noqa: F401
except ImportError:
    raise SystemExit(
        "ERROR: este script requiere el entorno pga-core.\n"
        "Ejecuta con:\n"
        "  ~/.pga/envs/pga-core/bin/python " + __file__
    )

sys.path.insert(0, str(Path(__file__).parent))
from genotype_simulator import simulate_snv_genotypes, simulate_sv_genotypes
from pedigree_model import build_standard_pedigree
from sv_catalog import default_catalog
from publicgenomicagent.knowledge.gnomad import fetch_region, mid_freq

ROI_CHROM = "2"
ROI_START = 110_000_000
ROI_END = 110_025_000
ROI_LEN = ROI_END - ROI_START
SEED = 42
CONTIG = f"chr{ROI_CHROM}_roi"

BCFTOOLS = str(Path("~/.pga/envs/pga-hts/bin/bcftools").expanduser())


def to_local(genomic_pos: int) -> int:
    """Convierte posición genómica 1-based a posición local 1-based del FASTA."""
    return genomic_pos - ROI_START


def load_mid_variants():
    raw = fetch_region(ROI_CHROM, ROI_START, ROI_END)
    out = []
    for v in raw:
        m = mid_freq(v)
        if m and m[0] >= 1:
            v["_af_mid"] = m[0] / m[1]
            out.append(v)
    return out


def write_plain(path, pedigree, snv_genos, sv_genos, snv_vars, sv_cat):
    samples = pedigree.order_for_vcf()
    records = []

    for v in snv_vars:
        pos_local = to_local(v["pos"])
        gts = [f"{snv_genos[i][v['variant_id']][0]}/{snv_genos[i][v['variant_id']][1]}" for i in samples]
        af = v.get("_af_mid", 0.0)
        info = f"AF_MID={af:.6f}"
        line = "\t".join([
            CONTIG, str(pos_local), v["variant_id"],
            v["ref"], v["alt"], ".", "PASS", info, "GT",
        ] + gts)
        records.append((pos_local, line))

    for sv in sv_cat:
        pos_local = to_local(sv.pos)
        end_local = to_local(sv.end)
        gts = [f"{sv_genos[i][sv.id][0]}/{sv_genos[i][sv.id][1]}" for i in samples]
        info = f"SVTYPE={sv.svtype};END={end_local};SVLEN={sv.length}"
        line = "\t".join([
            CONTIG, str(pos_local), sv.id,
            "N", f"<{sv.svtype}>", ".", "PASS", info, "GT",
        ] + gts)
        records.append((pos_local, line))

    records.sort(key=lambda x: x[0])

    with path.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=PublicGenomicAgent-fixture\n")
        f.write(f"##contig=<ID={CONTIG},length={ROI_LEN}>\n")
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n')
        f.write('##ALT=<ID=DEL,Description="Deletion">\n')
        f.write('##ALT=<ID=INS,Description="Insertion">\n')
        f.write('##ALT=<ID=DUP,Description="Duplication">\n')
        f.write('##ALT=<ID=INV,Description="Inversion">\n')
        f.write('##INFO=<ID=AF_MID,Number=1,Type=Float,Description="Allele frequency in gnomAD Middle Eastern population">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        f.write("\t".join(header_cols + samples) + "\n")
        for _, line in records:
            f.write(line + "\n")


def compress_and_index(plain, out_vcf):
    subprocess.run([BCFTOOLS, "view", str(plain), "-Oz", "-o", str(out_vcf)], check=True)
    subprocess.run([BCFTOOLS, "index", "-t", str(out_vcf)], check=True)


def write_fam(path, pedigree):
    with path.open("w") as f:
        for iid in pedigree.order_for_vcf():
            ind = pedigree.individuals[iid]
            pid = ind.father or "0"
            mid = ind.mother or "0"
            ph = "2" if ind.affected else "1"
            f.write(f"{iid}\t{iid}\t{pid}\t{mid}\t{ind.sex}\t{ph}\n")


def write_expected(path, snv_vars, sv_cat, sv_genos):
    payload = {
        "roi": {"chrom": ROI_CHROM, "start": ROI_START, "end": ROI_END},
        "contig": CONTIG,
        "snv_count": len(snv_vars),
        "sv_catalog": [
            {"id": s.id, "svtype": s.svtype,
             "pos_genomic": s.pos, "pos_local": to_local(s.pos),
             "length": s.length}
            for s in sv_cat
        ],
        "expected_findings": {
            "SV002_de_novo_C2": {
                "C2": list(sv_genos["C2"]["SV002"]),
                "F1": list(sv_genos["F1"]["SV002"]),
                "M1": list(sv_genos["M1"]["SV002"]),
            },
            "SV005_recessive_C3": {
                "C3": list(sv_genos["C3"]["SV005"]),
                "F1": list(sv_genos["F1"]["SV005"]),
                "M1": list(sv_genos["M1"]["SV005"]),
            },
        },
    }
    path.write_text(json.dumps(payload, indent=2))


def main():
    if len(sys.argv) != 2:
        print("Uso: make_cohort_fixture.py <outdir>", file=sys.stderr)
        sys.exit(2)
    outdir = Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(SEED)
    ped = build_standard_pedigree()
    snvs = load_mid_variants()
    cat = default_catalog(ROI_CHROM, ROI_START)
    snv_g = simulate_snv_genotypes(ped, snvs, rng)
    sv_g = simulate_sv_genotypes(ped, cat, rng)

    plain = outdir / "cohort.vcf"
    vcf = outdir / "cohort.vcf.gz"
    write_plain(plain, ped, snv_g, sv_g, snvs, cat)
    compress_and_index(plain, vcf)
    plain.unlink()
    write_fam(outdir / "pedigree.fam", ped)
    write_expected(outdir / "expected.json", snvs, cat, sv_g)

    print(f"OK: {vcf}")
    print(f"SNVs: {len(snvs)}  SVs: {len(cat)}  samples: {len(ped.individuals)}")
    print(f"Contig local: {CONTIG} (posiciones 1..{ROI_LEN})")


if __name__ == "__main__":
    main()
