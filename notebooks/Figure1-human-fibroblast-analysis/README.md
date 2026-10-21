# Figure 1 — human dermal fibroblasts, young (10 yo) vs old (75 yo)

Age-specific interchromosomal contacts (LAS) and pathway colocalization in human dermal fibroblasts, with expression in
the Hi-C donors and a 133-donor population cohort (Fleischer et al., GSE113957).

## Notebooks (run in this order, from this folder)
1. `fibroblast_cohort.ipynb` (needs pydeseq2 + anndata) — the population cohort → `data/`: DESeq2 of
   older (61–96 yo) vs younger (1–13 yo) healthy donors and GSEA (`cohort_60plus_vs_young_{deseq,gsea}.csv`), and the
   per-donor interferon-γ score (`cohort_isg_scores.csv`). ~1 min.
2. `fibroblast_analysis.ipynb` — every other computation → `data/`: pooled Hi-C maps (`fig1a_pooled_maps.npz`),
   pathway colocalization (`fig1_pathway_colocalization_tss.csv`; shared `../coloc_stats.py`), locus-matched permutation
   test (`locus_matched_fibroblast.csv`; shared `../locus_null.py`), donor GSEA (`75_vs_10_yo_gseapy/`; shared `../gsea.py`),
   the cohort IFN-γ score–age correlation (`cohort_ifng_age_spearman.csv`) and k-cores (shared `../kcore.py`;
   `kcores/<core>_{genes,states}.csv`, `kcores/kcore_summary.csv`). ~1–2 min.
3. `fibroblast_figures.ipynb` — reads only `data/` and draws every panel → `figures/<panel>.{svg,png}`, one
   `###` section per final panel. Seconds.
4. `fibroblast_pages.ipynb` — assembles the pages, cropped tight → `figures/Figure1.{pdf,png,svg}`,
   `figures/Figure1S.*`. Rule: rows fill the page width; main text prints at ≥ 7 pt where the layout allows (Figure 1S
   ~6.6 pt).

`figstyle.py` — Arial, editable SVG text. To run the notebooks headless, see the `nbclient` example in
`../Figure4-muscle/README.md`.

## Panel → notebook section → file
Each panel has a `###` section in `fibroblast_figures.ipynb` under the `##` section below; the computations are under the
same `##` sections, in the same order, in `fibroblast_analysis.ipynb` (cohort expression and scores in
`fibroblast_cohort.ipynb`).

| Panel | Content | `##` section | Panel file (`figures/`) |
|---|---|---|---|
| Figure 1a | Hi-C → LAS → age-specific contacts | Hi-C processing | `fig1a_hic_processing` |
| Figure 1b | colocalized pathways (hypergeometric FDR < 0.05) | Pathway colocalization | `fig1b_colocalization` |
| Figure 1c | colocalization vs expression, Hi-C donors | Colocalization vs expression | `fig1c_concordance_donors` |
| Figure 1d | population-cohort IFN-γ score vs age | Population cohort interferon score | `fig1d_cohort_ifng_age` |
| Figure 1e | IFN-γ old-specific 3-core | State-specific k-cores | `fig1e_ifng_core` |
| Supplementary 1Sa | locus-matched permutation test | Locus-matched permutation test | `fig1S_a_locus_matched` |
| Supplementary 1Sb | donor GSEA (75 vs 10 yo) | Pathway expression (GSEA) | `fig1S_b_expression_donors` |
| Supplementary 1Sc | donor colocalization vs cohort expression | Colocalization vs expression | `fig1S_c_concordance_cohort` |
| Supplementary 1Sd | cohort GSEA (61–96 vs 1–13 yo) | Pathway expression (GSEA) | `fig1S_d_expression_cohort` |
| Supplementary 1Se, 1Sf | young-specific cores (senescence inhibitors, ECM) | State-specific k-cores | `fig1S_e_seninh_core`, `fig1S_f_ecm_core` |

Gene assignment: one 250 kb bin per protein-coding gene (its TSS; `../../metadata/hg38/bins_protein_coding_genes_tss_hg38_250000.csv`).
Contact-state colours as Figures 4–6: old red, young blue, shared grey `#B0B0B0`, none `#F2F2F2`.

## Primary inputs (stay under `/data/land/3d-inter-paper-data/`)
- Age-specific LAS graph `networkx_graphs/G_all.pkl`; normalized / LAS Hi-C coolers `normalized_hic_cooler/`.
- Donor bulk RNA-seq counts `bulk-rna-data/RNAseq/old_vs_young_bulk.csv`; donor DESeq2 results
  `bulk-rna-data/old_vs_young_deseq_annotated.csv` (donor GSEA).
- Population cohort counts `jana-aging/rna_counts/pyximport_processed.h5ad` (133 healthy donors + 10 HGPS, excluded) and
  Jana's deposited Group4 vs Group1 DESeq2 result `jana-aging/deseq-young-vs-old-results/` (processing check).

## Notes
- The 60-pathway gene sets are read from `../../metadata/hallmark_and_custom_pathways_unfiltered_aug_2025.{json,gmt}` (the
  `.gmt` holds the same panel as gene-set lines, used by GSEA). `../../illustrator-figures/Figure1.ai` is the old
  (pre-2026-10) figure.
- Removed on 2026-10-20: the legacy Figure 1 / 1S / S2 notebooks, the global-enrichment test, the Group4/Group5 cohort GSEA,
  and the separate `gseapy.ipynb`, `cohort_60plus.py`, `cohort_isg_age.py` and `ifng_old_contact_genes.py` (folded into
  `fibroblast_analysis.ipynb` and `fibroblast_cohort.ipynb`; the exploratory IFN-α and old-contact-gene scores were dropped).
