"""Shared colocalization statistics for every dataset (fibroblast, muscle, PSC->CM, CD8).

The tests take unique interchromosomal gene-gene contact counts, split by state (state-A-only, state-B-only, shared),
for a gene set and for the genome-wide background; count_contacts computes them from a contact graph.

- count_contacts: unique gene-gene contacts per state, for the whole genome or restricted to a gene set.
- hypergeom_enrichment_pval: one-sided hypergeometric test that a set's state-X-only contacts exceed the genome rate
  (the test behind every bar panel / star).
- differential_fisher: two-sided Fisher's exact test of a set's B-only vs A-only contacts against the genome's
  B-only vs A-only contacts; the odds ratio equals fold_B / fold_A (reported as a complementary test).

Import from a notebook folder with:  sys.path.insert(0, "..");  import coloc_stats
"""
import itertools
from collections import defaultdict

import numpy as np
from scipy.stats import hypergeom, fisher_exact, false_discovery_control


def count_contacts(G, gene_attr, genes=None):
    """Unique gene-gene contacts per state ('edge_type' of each edge, e.g. state A / state B / 'Shared').
    G: networkx graph of interchromosomal contacts between 250 kb bins; each node lists its genes under `gene_attr`.
    genes: optional gene set; if given, only pairs of two of its genes count (the set's contacts). Every pair of
    genes on two contacting bins counts once per state, however many bin pairs carry it; a gene never pairs with itself."""
    node_genes = {n: list(gl) for n, gl in G.nodes(data=gene_attr) if gl}
    if genes is not None:
        keep = set(genes)
        node_genes = {n: [g for g in gl if g in keep] for n, gl in node_genes.items()}
        node_genes = {n: gl for n, gl in node_genes.items() if gl}
        G = G.subgraph(node_genes)                                  # only bins carrying genes of the set
    pairs = defaultdict(set)
    for b1, b2, state in G.edges(data="edge_type"):
        for g1, g2 in itertools.product(node_genes.get(b1, []), node_genes.get(b2, [])):
            if g1 != g2:
                pairs[state].add((g1, g2) if g1 < g2 else (g2, g1))
    return {state: len(p) for state, p in pairs.items()}


def hypergeom_enrichment_pval(k, M, n, N):
    """Hypergeometric P(X >= k): k state-specific contacts among a set's N, genome has n of M total contacts.
    Capped at 1 (sf + pmf can exceed 1 by floating-point rounding)."""
    return float(min(1.0, hypergeom.sf(k, M, n, N) + hypergeom.pmf(k, M, n, N)))


def differential_fisher(k_b, k_a, bg_b, bg_a):
    """Fisher's exact test (two-sided) of the set's B- vs A-specific contacts against the genome's.
    Returns (odds ratio = fold_B / fold_A, p); (nan, nan) if the set has no state-specific contacts."""
    if k_a + k_b == 0:
        return np.nan, np.nan
    odds_ratio, p = fisher_exact([[int(k_b), int(k_a)], [int(bg_b), int(bg_a)]])
    return float(odds_ratio), float(p)


def bh(pvals):
    """Benjamini-Hochberg FDR, ignoring NaNs."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    out[ok] = false_discovery_control(np.clip(p[ok], 0, 1), method="bh")
    return out

