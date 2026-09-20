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

# Variante de interés para demo de ganancia diagnóstica:
# AF_MID = 4.4%, heredada del padre (F1) en heterocigosis.
# C2 la porta en 0/1; F1 en 0/1; resto del pedigree 0/0.
FORCED_VARIANT_ID = "2-110008979-A-T"
FORCED_GENOTYPES = {
    FORCED_VARIANT_ID: {
        "F1": (0, 1),
        "C2": (0, 1),
    },
}
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



def _load_ref_seq_for_cluster():
    """Carga la referencia del ROI para que los SNVs sintéticos
    usen las bases reales de GRCh38 como REF."""
    src = Path.home() / ".pga/cache/reference/hg38_chr2_110000001_110025000.fa"
    if not src.exists():
        return None, 110_000_001
    lines = src.read_text().splitlines()
    seq = "".join(l for l in lines if not l.startswith(">"))
    return seq, 110_000_001


def add_mena_cluster(variants, center_pos, n=8, span=200, seed=7):
    """Añade n SNVs sintéticos MENA-específicos alrededor de center_pos.

    Simula un haplotipo MENA que GRCh38 no tiene. Usa las bases reales
    de GRCh38 como REF (leídas del FASTA cacheado).
    """
    import random as _rng
    rng = _rng.Random(seed)
    bases = "ACGT"

    ref_seq, ref_start = _load_ref_seq_for_cluster()
    if ref_seq is None:
        raise RuntimeError(
            "No se encuentra la referencia cacheada. "
            "Ejecuta primero el fetch de UCSC."
        )

    used = {v["pos"] for v in variants}
    new_variants = []
    attempts = 0
    while len(new_variants) < n and attempts < 500:
        attempts += 1
        offset = rng.randint(-span // 2, span // 2)
        if offset == 0:
            continue
        pos = center_pos + offset
        if pos in used:
            continue
        # Coordenada local en el FASTA cacheado
        local_idx = pos - ref_start
        if local_idx < 0 or local_idx >= len(ref_seq):
            continue
        ref = ref_seq[local_idx].upper()
        alt = rng.choice([b for b in bases if b != ref])
        vid = f"MENA_cluster_{pos}_{ref}_{alt}"
        new_variants.append({
            "variant_id": vid,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "_af_mid": rng.uniform(0.3, 0.6),
            "_is_synthetic_mena": True,
        })
        used.add(pos)

    return variants + new_variants

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

    # Añadir cluster MENA alrededor de la variante objetivo para forzar
    # reference bias. Las coordenadas de las variantes del cohort están
    # en GENÓMICAS (se convierten a locales en write_plain).
    # La variante objetivo es 2-110008979-A-T (genómica).
    target_pos_genomic = 110_008_979
    snvs = add_mena_cluster(snvs, center_pos=target_pos_genomic, n=8, span=200)

    cat = default_catalog(ROI_CHROM, ROI_START)

    # El cluster MENA debe estar presente en C2 y sus padres
    # (son SNVs comunes en MENA, no de novo)
    mena_cluster_ids = [v["variant_id"] for v in snvs if v.get("_is_synthetic_mena")]
    extra_forced = {}
    for vid in mena_cluster_ids:
        extra_forced[vid] = {
            "GF1": (0, 1), "GM1": (0, 1), "GF2": (0, 1), "GM2": (0, 1),
            "F1": (0, 1), "M1": (0, 1),
            "C1": (0, 1), "C2": (0, 1), "C3": (0, 1), "C4": (0, 1),
        }
    forced_all = {**FORCED_GENOTYPES, **extra_forced}

    snv_g = simulate_snv_genotypes(ped, snvs, rng, forced=forced_all)
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
