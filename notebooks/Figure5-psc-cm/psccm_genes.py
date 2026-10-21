"""Gene -> bin assignment for the PSC->CM LAS graphs: protein-coding genes by TSS (hg19, GENCODE v19), the same
convention as the fibroblast / muscle / CD8 analyses. The graphs' own `genes` attribute lists every gene overlapping
a bin; use_tss_genes() overwrites it so all downstream counting uses one TSS bin per protein-coding gene."""
import pandas as pd
import networkx as nx

TSS_BINS = "../../metadata/hg19/bins_df_protein_coding_genes_with_tss_hg19_250000.csv"   # relative to this notebook folder
_tss_genes = None


def tss_genes():
    """{node id: [protein-coding genes whose TSS falls in that 250 kb bin]} (built by metadata/hg19/build_hg19_tss_annotation.py)."""
    global _tss_genes
    if _tss_genes is None:
        col = pd.read_csv(TSS_BINS, index_col=0)["protein_coding_genes_with_tss"]
        _tss_genes = {i: s.split(",") if isinstance(s, str) and s else [] for i, s in col.items()}
    return _tss_genes


def use_tss_genes(G):
    """Overwrite each node's `genes` attribute with its TSS-assigned protein-coding genes; returns G."""
    genes = tss_genes()
    nx.set_node_attributes(G, {n: genes.get(n, []) for n in G.nodes}, "genes")
    return G
