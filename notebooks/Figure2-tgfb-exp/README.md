# Figure 2 — TGF-β stimulation of young and old fibroblasts (RNA-seq + ATAC-seq)

Young (10 yo) and old (75 yo) fibroblasts under mechanical tension, control or TGF-β, two replicates each (bulk RNA-seq,
GSE300230; ATAC-seq, GSE300231). Asks whether genes that colocalize in young fibroblasts (Figure 1) respond to TGF-β
differently from other genes of the same pathway, in young vs old cells.

## Notebooks (run in this order, from this folder)
1. `tgfb_analysis.ipynb` — every computation → `data/`: expression and DE table (`expression_tgfb.csv`), contact classes
   per pathway (`pathway_gene_contacts.csv`, `pathway_subsets.csv`), Figure 2c Mann–Whitney tests
   (`fig2c_mannwhitney.csv`), gene-level contact graphs (`gene_graphs/<pathway>_<subset>_{nodes,edges}.csv`), promoter
   annotation (`promoters_tss_1kb.bed`, `consensus_peak_annotated.bed`; needs `bedtools`),
   promoter accessibility (`fig2d_promoter_accessibility.csv`), promoter-peak changes (`fig2Sa_c_promoter_peak_de*.csv`) and
   the gene-density control (`fig2Sb_gene_groups.csv`, `fig2Sb_mannwhitney.csv`). ~30 s.
2. `tgfb_figures.ipynb` — reads only `data/` and draws every panel → `figures/<panel>.svg` (2d also `.pdf`, 2Sb also
   `.png`), one `###` section per final panel (needs `pygraphviz` and `pdftocairo`). ~15 s.
3. `tgfb_pages.ipynb` — assembles `Figure2`, `Figure2Sa` and `Figure2Sb`, cropped tight → PDF/PNG in `figures/`, page SVGs
   in `../../illustrator-figures/`.

To run the notebooks headless, see the `nbclient` example in `../Figure4-muscle/README.md`.

## Panel → notebook section → file
Each panel has a `###` section in `tgfb_figures.ipynb` under the `##` section below; the computations are under the same
`##` sections, in the same order, in `tgfb_analysis.ipynb`.

| Panel | Content | `##` section | Panel file (`figures/`) |
|---|---|---|---|
| Figure 2a | experiment schematic (hand-drawn vector) | — | `data/fig2a_schematic.svg` |
| Figure 2b | ECM pathway contact graphs (young-only / other / old-only / other) | Gene-level contact graphs | `fig2b_ecm_contact_graphs` |
| Figure 2c | TGF-β response by contact class, young and old cells | TGF-β response by contact class | `fig2c_tgfb_response` |
| Figure 2d | promoter accessibility, young-specific vs other genes | Promoter accessibility | `fig2d_promoter_accessibility` |
| Supplementary 2Sa, panel a | expression histograms per condition | Expression | `fig2Sa_a_expression_histograms` |
| Supplementary 2Sa, panel b | contact graphs of TGF-β signaling, senescence inhibitors, MYC targets | Gene-level contact graphs | `fig2Sa_b_contact_graphs` |
| Supplementary 2Sa, panel c | promoter-peak change under TGF-β | Promoter peak changes under TGF-β | `fig2Sa_c_promoter_peak_de` |
| Supplementary 2Sb, panel a | locus gene density by contact group | Locus gene density and contact groups | `fig2Sb_a_gene_density` |
| Supplementary 2Sb, panel b | TGF-β response by contact group, young vs old | Locus gene density and contact groups | `fig2Sb_b_response_by_contact_group` |

`fig2_legend_edge_types.svg` is the edge-type legend placed next to 2b and 2Sa b by the pages notebook.

Definitions: a gene is active if its mean count over the two replicates exceeds 100 in any of the four tension conditions
(genes on chrX/Y/M excluded, as the Hi-C analysis). Each gene sits on the 250 kb bin holding its TSS
(`../../metadata/hg38/bins_protein_coding_genes_tss_hg38_250000.csv`; same rule as Figures 1 and 4–6). A gene forms
young-specific (old-specific) contacts if its bin has a young-only (old-only) contact with the bin of another active gene
of the same pathway. Promoters = TSS ± 1 kb.

## Inputs
- In `data/`: RNA-seq counts `tgfb-force-exp-counts.csv` and DESeq2 results `deseq/<cond1>_vs_<cond2>_deseq_results.csv`;
  consensus ATAC peaks of the 8 samples `consensus_peak-8samples.bed` and per-sample MACS2 peaks `atac-seq/raw/`.
- Under `/data/land/3d-inter-paper-data/`: the age-specific contact graph `networkx_graphs/G_all.pkl`, and the DESeq2
  results of the ATAC peaks `yawen-atac-data/de-peaks/{YRT_vs_YRC,ORT_vs_ORC}_deseq_peaks_results.csv`.
- Shared: `../../metadata/hallmark_and_custom_pathways_unfiltered_aug_2025.json`, `../../metadata/hg38/`.
- Raw data are in GEO: GSE300230 (RNA-seq raw counts) and GSE300231 (ATAC-seq MACS2 peaks = `data/atac-seq/raw/`).

## Notes
- The gene-graph layouts (graphviz) follow the order in which genes enter the graph, which depends on Python's per-session
  set ordering, so 2b and 2Sa b are laid out differently on each fresh run of `tgfb_analysis.ipynb` (same nodes, edges and
  colours). Run with `PYTHONHASHSEED=0` to reproduce a layout exactly; the saved node/edge order in `data/gene_graphs/`
  is what the figures notebook draws.
- Cleaned up 2026-10-20: `tgfb_figures.ipynb` (computations + plots) was split into analysis / figures; the old-pathway
  (IL2-STAT5, bile acid) panels, the two-column graph version, the deprecated ATAC peak-count violins and unused legends
  were dropped; `analysis.ipynb` (exploratory sweeps; its promoter recipe now lives in `tgfb_analysis.ipynb`),
  `atac_preprocessing.ipynb` and their unused outputs were removed. Supplementary Figure 2Sb (gene-density control) added; the first supplementary page was renamed 2S → 2Sa
  (as Figures 4 and 5), and its letter "b." is no longer covered by the edge-type legend.
  2Sa a is now drawn by the notebook (its save line had been commented out); it shows the same data as the previous
  file, which was saved with an older matplotlib (3.7.3), so the bars differ by ~0.1% in height.
