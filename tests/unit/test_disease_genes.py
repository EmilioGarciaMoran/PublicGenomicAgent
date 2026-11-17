"""Fenotipo × genotipo: OMIM -> genes (ficheros de LIRICAL) y cruce con candidatas."""
from __future__ import annotations

from pathlib import Path

from publicgenomicagent.knowledge.disease_genes import concordance, load_disease_genes, read_lirical_tsv

MINI = Path(__file__).parents[1] / "fixtures" / "lirical_mini"


def test_omim_to_genes_and_aliases():
    dg = load_disease_genes(MINI)
    assert dg.omim_to_genes["OMIM:256100"] == {"NPHP1"}           # filas phenotype con GeneID
    assert "OMIM:607100" not in dg.omim_to_genes                  # filas gene: no son enfermedades
    assert "OMIM:999999" not in dg.omim_to_genes                  # sin GeneID
    assert dg.symbol("slsn1") == "NPHP1" and dg.symbol("MKS3") == "TMEM67"
    assert dg.symbol("GENSIM") is None


def test_concordance_side_by_side_with_ties():
    dg = load_disease_genes(MINI)
    ranking = read_lirical_tsv(MINI / "ranking.tsv")
    cands = [{"gene": "OLDNPHP", "model": "recessive homozygous", "variant": "chr2 110008979 A T"},
             {"gene": "GENSIM", "model": "recessive homozygous", "variant": "chr1 500 C G"},
             {"gene": "GENEX", "model": "compound het (trans)", "variant": "a + b"}]
    res = concordance(cands, ranking, dg, top_n=5)
    first = res["candidates"][0]
    assert first["gene_symbol"] == "NPHP1" and first["gene_status"] == "ranked"
    b = first["best_disease"]
    assert (b["rank"], b["disease_curie"], b["tied_with"]) == (4, "OMIM:256100", 3)   # 4 con LR 5.253
    status = {c["gene"]: c["gene_status"] for c in res["candidates"]}
    assert status == {"OLDNPHP": "ranked", "GENSIM": "not_a_gene", "GENEX": "no_ranked_disease"}
    without = {t["disease_curie"]: t["genes"] for t in res["top_diseases_without_candidate"]}
    assert res["n_ranked_diseases"] == 8525          # de pretestprob 1/8525, no las 6 filas del TSV
    assert without == {"OMIM:606996": ["NPHP4"], "OMIM:606966": ["NPHP4"], "OMIM:619113": ["TMEM67"]}
