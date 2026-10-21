"""Shared locus-matched permutation test for state-specific pathway colocalization (every dataset).

Complements coloc_stats.hypergeom_enrichment_pval, which treats a set's gene pairs as independent draws. Here a set's
fold enrichment of state-specific contacts is compared with random sets of look-alike loci:

- each gene sits on the 250 kb bin holding its TSS; a set becomes a list of loci weighted by its genes per locus, and
  its statistic is the share of its contacting interchromosomal gene pairs that are state-X-only (pairs between loci
  i, j count w_i * w_j; fold = share / genome-wide share, identical to the hypergeometric panels);
- null: every locus is replaced by a random gene-bearing locus from the same stratum, drawn without replacement with
  probability proportional to its genes, weights kept. Strata = chromosome x gene density (1, 2-3, 4+ protein-coding
  genes) x expression in state X (0 = every gene silent, else quartile of the bin's maximum gene expression among
  expressed bins). Without RNA (expr=None) the expression term is dropped. A stratum smaller than the number of loci
  to draw falls back to the whole chromosome;
- p = (null shares >= observed + 1) / (n_draws + 1), one-sided per state; BH per state over the gene-set family.

Deterministic: each set uses seed 1000 + its position in the family, and strata are visited in sorted order.
Usage (from a notebook folder):
    sys.path.insert(0, ".."); import locus_null
    locus_null.prepare(G, bins, states=("young", "old"), expr={"young": s1, "old": s2}, gene_col="protein_coding_genes_with_tss")
    table = locus_null.test_family(gene_sets)          # DataFrame: set, n_genes, n_loci, pairs, fold_*, locus_p_*, locus_fdr_*
"""
from multiprocessing import get_context

import numpy as np
import pandas as pd

import coloc_stats

DENSITY = [1, 2, 4]          # gene-density classes: 1, 2-3, 4+ protein-coding TSS genes per locus
CHUNK = 250                  # null draws scored per vectorized batch
CTX = {}                     # dataset arrays, set by prepare() before the worker pool forks


def prepare(G, bins, states, expr=None, gene_col="genes"):
    """Build the contact-state matrix, locus strata and genome-wide background for one graph into CTX.
    G: networkx graph whose nodes are row positions of `bins` and whose edges carry 'edge_type' (state name or 'Shared').
    bins: DataFrame (one row per bin, in node order) with 'chrom' and `gene_col` (comma-separated genes).
    states: (state A, state B) as in the edge types. expr: None, or {state: Series of expression by gene}."""
    n = len(bins)
    a, b = states
    code = {a: 1, b: 2, "Shared": 3}
    S = np.zeros((n, n), np.int8)
    for u, v, t in G.edges(data="edge_type"):
        S[u, v] = S[v, u] = code[t]
    chrom = bins["chrom"].astype(str).to_numpy()
    genes_col = bins[gene_col].to_numpy()
    bin_genes = {i: v.split(",") for i, v in enumerate(genes_col) if isinstance(v, str) and v}
    n_genes = np.zeros(n)
    for i, gs in bin_genes.items():
        n_genes[i] = len(gs)
    loci = np.array(sorted(bin_genes))
    strata, pools = {}, {}
    for st in (a, b):
        if expr is not None:                                       # bin max expression; quartile cuts over expressed bins
            bmax = {i: expr[st].reindex(bin_genes[i]).fillna(0).max() for i in loci}
            cuts = np.quantile([v for v in bmax.values() if v > 0], [.25, .5, .75])

        def key(i):
            k = [chrom[i], int(np.searchsorted(DENSITY, n_genes[i], side="right"))]
            if expr is not None:
                k.append(0 if bmax[i] == 0 else 1 + int(np.searchsorted(cuts, bmax[i], side="right")))
            return tuple(k)
        strata[st] = {i: key(i) for i in loci}
        pools[st] = {}
        for i, k in strata[st].items():
            pools[st].setdefault(k, []).append(i)
        pools[st] = {k: np.array(v) for k, v in pools[st].items()}
    iu = np.triu_indices(n, 1)
    wp, sv = n_genes[iu[0]] * n_genes[iu[1]], S[iu]
    bg = {s: float(wp[sv == c].sum()) for s, c in code.items()}
    CTX.clear()
    CTX.update(S=S, chrom=chrom, n_genes=n_genes, loci=loci, strata=strata, pools=pools, code=code, bg=bg,
               bg_total=sum(bg.values()), gene_bin={g: i for i, gs in bin_genes.items() for g in gs}, states=(a, b))


def set_loci(genes):
    """Gene set -> (loci, genes per locus)."""
    s = pd.Series([CTX["gene_bin"][g] for g in set(genes) if g in CTX["gene_bin"]]).value_counts().sort_index()
    return s.index.to_numpy(), s.to_numpy().astype(float)                  # loci in bin order (deterministic)


def shares(idx, w, code):
    """State share and number of contacting interchromosomal gene pairs, for each row of loci idx (draws x loci)."""
    S, chrom = CTX["S"], CTX["chrom"]
    sub = S[idx[:, :, None], idx[:, None, :]]
    W = np.outer(w, w)
    W[chrom[idx[0]][:, None] == chrom[idx[0]][None, :]] = 0        # interchromosomal pairs only
    tot = ((sub > 0) * W).sum(axis=(1, 2)) / 2
    return np.where(tot > 0, ((sub == code) * W).sum(axis=(1, 2)) / 2 / np.maximum(tot, 1), 0.0), tot


def null_shares(loci, w, state, n_draws, rng):
    """Null state shares: loci replaced within their strata (probability ~ genes per locus), weights kept."""
    keys = [CTX["strata"][state][i] for i in loci]
    out = []
    for start in range(0, n_draws, CHUNK):
        n = min(CHUNK, n_draws - start)
        idx = np.empty((n, len(loci)), dtype=np.int64)
        for key in sorted(set(keys)):                              # sorted: identical random stream on every run
            pos = [j for j, k in enumerate(keys) if k == key]
            pool = CTX["pools"][state][key]
            if len(pool) < len(pos):                               # tiny stratum: fall back to the chromosome
                pool = CTX["loci"][CTX["chrom"][CTX["loci"]] == key[0]]
            # weighted sampling without replacement (exponential keys), probability ~ genes per locus
            draws = np.argsort(-np.log(rng.random((n, len(pool)))) / CTX["n_genes"][pool], axis=1)[:, :len(pos)]
            idx[:, pos] = pool[draws]
        out.append(shares(idx, w, CTX["code"][state])[0])
    return np.concatenate(out)


def test_set(args):
    """Fold enrichment, hypergeometric p and locus-matched p for one gene set, both states."""
    name, genes, n_draws, seed = args
    rng = np.random.default_rng(seed)
    loci, w = set_loci(genes)
    row = {"set": name, "n_genes": int(w.sum()), "n_loci": len(loci)}
    if len(loci) < 3:
        return row
    for st in CTX["states"]:
        obs, tot = shares(loci[None, :], w, CTX["code"][st])
        obs, tot = float(obs[0]), float(tot[0])
        row["pairs"] = tot
        row[f"fold_{st}"] = obs / (CTX["bg"][st] / CTX["bg_total"]) if tot else np.nan
        row[f"hyper_p_{st}"] = (coloc_stats.hypergeom_enrichment_pval(round(obs * tot), CTX["bg_total"], CTX["bg"][st], tot)
                                if tot else np.nan)
        null = null_shares(loci, w, st, n_draws, rng)
        row[f"locus_p_{st}"] = ((null >= obs).sum() + 1) / (n_draws + 1) if tot else np.nan
    return row


def test_family(gene_sets, n_draws=10_000, processes=16):
    """Every set of a family (dict name -> genes; seed = 1000 + position), with BH FDR per state. Call prepare() first."""
    jobs = [(name, genes, n_draws, 1000 + i) for i, (name, genes) in enumerate(gene_sets.items())]
    with get_context("fork").Pool(processes) as pool:              # workers inherit CTX
        res = pd.DataFrame(pool.map(test_set, jobs))
    for st in CTX["states"]:
        for t in ("hyper", "locus"):
            res[f"{t}_fdr_{st}"] = coloc_stats.bh(res[f"{t}_p_{st}"])
    return res
