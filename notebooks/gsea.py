"""Shared preranked GSEA for every RNA-seq comparison (fibroblast, muscle, PSC->CM, CD8).

- deseq_rank: genes ranked by the DESeq2 Wald statistic (signed), the standard input for preranked GSEA.
- prerank: GSEApy preranked GSEA with the paper's fixed settings (1000 permutations, seed 6, sets of 5-5000 genes).

Import from a notebook folder with:  sys.path.insert(0, "..");  import gsea
"""
import pandas as pd
import gseapy as gp


def deseq_rank(de, genes=None):
    """DESeq2 results -> ranking (gene, score = Wald stat), highest first; `genes` = gene-name column (default: index).
    Genes without a statistic are dropped; for a duplicated gene name the copy with the largest |stat| is kept."""
    names = de.index if genes is None else de[genes]
    r = pd.DataFrame({"gene": names.values, "score": de["stat"].values}).dropna()
    r = r.loc[r.score.abs().sort_values(ascending=False).index].drop_duplicates("gene")
    return r.sort_values("score", ascending=False).reset_index(drop=True)


def prerank(rnk, gene_sets):
    """Preranked GSEA -> GSEApy's results table (Term, ES, NES, NOM p-val, FDR q-val, ...) with numeric statistics."""
    res = gp.prerank(rnk=rnk, gene_sets=gene_sets, min_size=5, max_size=5000, permutation_num=1000, seed=6,
                     threads=1, outdir=None, no_plot=True).res2d
    for col in ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]:
        res[col] = pd.to_numeric(res[col])
    return res
