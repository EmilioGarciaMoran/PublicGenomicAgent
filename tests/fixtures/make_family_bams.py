"""BAMs sintéticos de una familia para probar la verificación familiar.

Los haplotipos se transmiten con recombinación (ver `gamete`).

Genera una referencia aleatoria, un panel de SNV comunes y, por cada
individuo, un BAM con lecturas a profundidad fija cuyos alelos siguen
haplotipos transmitidos de forma mendeliana. Permite "intercambiar" un
BAM: escribir con el SM de un individuo las lecturas de otro.
Requiere pysam.
"""
from __future__ import annotations

import random
from pathlib import Path

SWAP = {"A": "G", "C": "T", "G": "A", "T": "C"}


def make_family(
    out: Path,
    parents: dict[str, tuple[str | None, str | None]],
    *,
    n_sites: int = 300,
    depth: int = 25,
    seed: int = 11,
    bam_as: dict[str, str] | None = None,
    recomb: float = 0.02,
) -> dict:
    """parents: sample -> (padre, madre). bam_as: sample -> individuo cuyas lecturas lleva.

    Devuelve {"ref": Path, "sites": Path, "bams": {sample: Path}}.
    """
    import pysam

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    length = 500 + 200 * n_sites + 10_000
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    ref = out / "ref.fa"
    ref.write_text(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, length, 60)) + "\n")
    pysam.faidx(str(ref))

    sites = [(p, seq[p - 1], SWAP[seq[p - 1]], rng.uniform(0.2, 0.5))
             for p in range(500, 500 + 200 * n_sites, 200)]
    sites_tsv = out / "sites.tsv"
    sites_tsv.write_text("".join(f"chr1\t{p}\t{r}\t{a}\n" for p, r, a, _ in sites))

    haps: dict[str, list[list[int]]] = {}

    def get(s: str) -> list[list[int]]:
        if s not in haps:
            f, m = parents[s]
            if f is None:
                haps[s] = [[int(rng.random() < af) for *_, af in sites] for _ in range(2)]
            else:
                haps[s] = [gamete(get(f)), gamete(get(m))]
        return haps[s]

    def gamete(parent: list[list[int]]) -> list[int]:
        """Haplotipo transmitido con recombinación (prob. `recomb` de cambio por sitio).

        Sin recombinación, dos hermanos podrían heredar haplotipos opuestos en
        todos los sitios y parecer no emparentados, algo imposible en una
        familia real a escala del genoma.
        """
        h = rng.randrange(2)
        out = []
        for k in range(len(sites)):
            if k and rng.random() < recomb:
                h = 1 - h
            out.append(parent[h][k])
        return out

    bams: dict[str, Path] = {}
    for sample in parents:
        source = (bam_as or {}).get(sample, sample)
        h = get(source)
        hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
               "SQ": [{"SN": "chr1", "LN": length}],
               "RG": [{"ID": "rg", "SM": sample}]}
        reads = []
        n = 0
        for k, (p, _r, a, _af) in enumerate(sites):
            for _ in range(depth):
                start = p - 1 - rng.randint(0, 99)
                s = list(seq[start:start + 100])
                if h[rng.randrange(2)][k]:
                    s[p - 1 - start] = a
                x = pysam.AlignedSegment()
                x.query_name = f"{sample}_{n}"
                n += 1
                x.query_sequence = "".join(s)
                x.flag = 0
                x.reference_id = 0
                x.reference_start = start
                x.mapping_quality = 60
                x.cigar = [(0, 100)]
                x.query_qualities = pysam.qualitystring_to_array("I" * 100)
                x.set_tag("RG", "rg")
                reads.append(x)
        reads.sort(key=lambda r: r.reference_start)
        path = out / f"{sample}.bam"
        with pysam.AlignmentFile(str(path), "wb", header=hdr) as o:
            for x in reads:
                o.write(x)
        pysam.index(str(path))
        bams[sample] = path
    return {"ref": ref, "sites": sites_tsv, "bams": bams}
