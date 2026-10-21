"""Shared state-specific k-core analysis for the module heatmaps (Figures 1, 4, 5 and 6).

Among a gene set, two genes are linked when their interchromosomal contact is specific to one state (e.g. old-only);
the k-core of that graph (default: the densest) is the module, shown as a gene x gene contact-state heatmap with
chromosome strips. Each gene sits on the 250 kb bin holding its TSS. A gene pair's state is the majority edge type over
the genes' bins (tie -> the shared label); with one TSS bin per gene it is simply the edge type between the two bins.
Core genes are returned in genomic order (chromosome, then bin).

Usage (from a notebook folder):
    sys.path.insert(0, ".."); import kcore
    edge_type = nx.get_edge_attributes(G, "edge_type")                  # {(bin, bin): state}
    genes, k = kcore.state_kcore(edge_type, gene_bins, pathway_genes, "Old", bin_chrom)     # densest core
    genes, k = kcore.state_kcore(edge_type, gene_bins, pathway_genes, "Old", bin_chrom, k=3)  # a fixed k
    mat = kcore.state_matrix(edge_type, gene_bins, genes, {"Young": 1, "Old": 2, "Shared": 3})
    kcore.chrom_strips(ax, genes, [bin_chrom[min(gene_bins[g])] for g in genes], fs=10)
gene_bins: {gene: set of its bins}; bin_chrom: {bin: chromosome number (1-22)}.
"""
import itertools
from collections import Counter

import networkx as nx
import numpy as np
from matplotlib.colors import to_rgb
from mpl_toolkits.axes_grid1 import make_axes_locatable

# paper-wide chromosome palette, chr1-22
CHROM_COLORS = dict(zip(range(1, 23), [
    "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928", "#1B9E77", "#E7298A", "#66A61E", "#E6AB02",
    "#666666", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#83d5af", "#ff69b4", "#b249d5", "#7edc45",
    "#c497cf", "#594e50"]))


def pair_state(edge_type, gene_bins, g1, g2, shared="Shared"):
    """Contact state of a gene pair: the majority edge type over the genes' bins (tie -> `shared`), else None."""
    states = []
    for a in gene_bins[g1]:
        for b in gene_bins[g2]:
            if a == b:
                continue
            t = edge_type.get((a, b)) or edge_type.get((b, a))
            if t:
                states.append(t)
    counts = Counter(states).most_common()
    if not counts:
        return None
    if len(counts) > 1 and counts[0][1] == counts[1][1]:
        return shared
    return counts[0][0]


def state_kcore(edge_type, gene_bins, genes, state, bin_chrom, k=None, shared="Shared"):
    """k-core of the `state`-specific colocalization graph among the genes (default: the densest).
    Returns (core genes in genomic order, k); a fixed k keeps every gene with core number >= k."""
    pairs = [(g1, g2) for g1, g2 in itertools.combinations(genes, 2)
             if pair_state(edge_type, gene_bins, g1, g2, shared) == state]
    core = nx.core_number(nx.Graph(pairs))
    k = max(core.values()) if k is None else k
    module = [g for g, v in core.items() if v >= k]
    return sorted(module, key=lambda g: (bin_chrom[min(gene_bins[g])], min(gene_bins[g]))), k


def state_matrix(edge_type, gene_bins, genes, codes, shared="Shared"):
    """Gene x gene contact-state matrix: codes[state] per pair, 0 = no interchromosomal contact."""
    mat = np.zeros((len(genes), len(genes)), int)
    for i, j in itertools.combinations(range(len(genes)), 2):
        state = pair_state(edge_type, gene_bins, genes[i], genes[j], shared)
        mat[i, j] = mat[j, i] = codes.get(state, 0)
    return mat


def chrom_strips(ax, genes, chroms, fs, strip_size=0.22, min_number_fs=None):
    """Chromosome strips between the gene labels and a heatmap (left + bottom), numbered once per chromosome block.
    strip_size: strip width (inches); min_number_fs: lower bound for the block-number font (default fs - 2)."""
    div = make_axes_locatable(ax)
    left = div.append_axes("left", size=strip_size, pad=0.03)
    bottom = div.append_axes("bottom", size=strip_size, pad=0.03)
    number_fs = fs - 2 if min_number_fs is None else max(fs - 2, min_number_fs)
    rgb = np.array([to_rgb(CHROM_COLORS[c]) for c in chroms])
    left.imshow(rgb[:, None, :], aspect="auto")
    bottom.imshow(rgb[None, :, :], aspect="auto")
    starts = [0] + [i for i in range(1, len(chroms)) if chroms[i] != chroms[i - 1]] + [len(chroms)]
    for a, b in zip(starts[:-1], starts[1:]):
        r, g, bl = to_rgb(CHROM_COLORS[chroms[a]])
        tc = "white" if 0.299 * r + 0.587 * g + 0.114 * bl < 0.55 else "black"     # readable on dark strips
        mid = (a + b - 1) / 2                                                        # centre of this chromosome block
        for strip, xy in ((left, (0, mid)), (bottom, (mid, 0))):
            strip.text(*xy, str(chroms[a]), ha="center", va="center", fontsize=number_fs, color=tc)
    left.set_xticks([])
    left.set_yticks(range(len(genes)), genes, fontsize=fs)
    bottom.set_yticks([])
    bottom.set_xticks(range(len(genes)), genes, rotation=90, fontsize=fs)
    for a in (left, bottom):
        a.tick_params(length=0)
        for sp in a.spines.values():
            sp.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
