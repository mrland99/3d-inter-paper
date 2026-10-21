# Figure 5 — PSC → cardiomyocyte interchromosomal colocalization (GSE116862)

Human ESC (H9, day 0) vs ventricular cardiomyocytes (day 80), Zhang et al. 2019; 2 biological replicates per state.
Interchromosomal (*trans*) LAS pathway colocalization at 250 kb (hg19), paired with RNA-seq.

## Notebooks (run in this order, from this folder)
1. `psccm_expression.ipynb` (needs pydeseq2 0.4.12, gseapy 1.3.1) — recount3 counts → DESeq2 (day 80 vs day 0)
   → GSEA on the Wald statistic (shared `../gsea.py`) → `data/gsea_{hallmark,cardiac,stemness,celltype}.csv`. ~2 min.
2. `psccm_analysis.ipynb` — every other computation → `data/`: colocalization test of every gene set on the
   pooled graph and the 6 replicate-pair graphs (shared `../coloc_stats.py`), the joint FDR and gene-set tables, the
   per-state / Fisher supplementary tables, the CM-specific k-cores (shared `../kcore.py`) and the replicate agreement.
   ~4 min (measured 4–6 min).
3. `psccm_figures.ipynb` (needs the Arial font) — reads only `data/` and draws every panel →
   `figures/<panel>.{svg,png}`, one `###` section per final panel. ~10 s.
4. `psccm_pages.ipynb` (needs `cairosvg`) — assembles `Figure5`, `Figure5Sa` and `Figure5Sb`, cropped tight →
   editable SVG in `../../illustrator-figures/`, PDF/PNG in `figures/`.

Helpers: `psccm_genes.py` assigns each protein-coding gene to its TSS bin (`../../metadata/hg19/`, built by
`build_hg19_tss_annotation.py` there; same convention as the other figures); `figstyle.py` sets Arial and editable SVG text.

To run all notebooks headless, e.g. with `nbclient`:

```python
import nbformat
from nbclient import NotebookClient
KERNEL = "python3"                                       # any Jupyter kernel with the required packages
for nb_name in ["psccm_expression.ipynb", "psccm_analysis.ipynb", "psccm_figures.ipynb", "psccm_pages.ipynb"]:
    nb = nbformat.read(nb_name, as_version=4)
    NotebookClient(nb, kernel_name=KERNEL, timeout=3600,
                   resources={"metadata": {"path": "."}}).execute()   # run from this folder so data/ and figures/ resolve
    nbformat.write(nb, nb_name)
```

## Panel → notebook section → file
Each panel has a `###` section in `psccm_figures.ipynb` under the `##` section below; it is computed in the
`psccm_analysis.ipynb` section whose header names the panel (GSEA in `psccm_expression.ipynb`).

| Panel | Content | `##` section (figures) | Data (`data/`) | Panel file (`figures/`) |
|---|---|---|---|---|
| Figure 5a | differentiation schematic | Schematic | — | `fig5a_schematic` |
| Figure 5b | colocalization + expression of the 5 main cardiomyocyte sets | Colocalization and expression of the cardiomyocyte gene sets | `panel_table.csv` | `fig5b_coloc_expression` |
| Figure 5c | colocalization vs expression, 60-pathway panel + 5 cardiomyocyte sets | Colocalization vs expression | `colocalization.csv` + `supp_coloc_vs_expression_table.csv` | `fig5c_coloc_vs_expression` |
| Figure 5d | cell-type specificity, 178 PanglaoDB cell types | Colocalization vs expression | `colocalization.csv` + `gsea_celltype.csv` | `fig5d_celltype_specificity` |
| Figure 5e | Sarcomere organization (GO) CM-specific core cluster | CM-specific core clusters | `kcores/sarcomere_*.csv` | `fig5e_sarcomere_core_cluster` |
| Figure 5f | Cardiomyocytes (PanglaoDB) CM-specific core cluster | CM-specific core clusters | `kcores/cardiomyocyte_*.csv` | `fig5f_cardiomyocyte_core_cluster` |
| Supplementary 5Sa, panel a | replicate design | Robustness | — | `fig5Sa_a_replicate_design` |
| Supplementary 5Sa, panel b | cardiomyocyte sets in each replicate pair | Robustness | `colocalization.csv` (Spearman for the legend text: `replicate_agreement.csv`) | `fig5Sa_b_replicate_sets` |
| Supplementary 5Sa, panel c | all 14 published cardiomyocyte sets | Robustness | `supp_cardiac_sets_table.csv` | `fig5Sa_c_all_cardiac_sets` |
| Supplementary 5Sa, panel d | all 13 published stemness sets | Robustness | `supp_stemness_sets_table.csv` | `fig5Sa_d_all_stemness_sets` |
| Supplementary 5Sb, panel a | colocalized pathways of the 60-pathway panel | Pathway context | `supp_coloc_vs_expression_table.csv` | `fig5Sb_a_pathway_colocalization` |
| Supplementary 5Sb, panel b | differentially expressed pathways of the 60-pathway panel | Pathway context | `supp_coloc_vs_expression_table.csv` | `fig5Sb_b_pathway_expression` |

## Outputs in `data/`
- `gsea_{hallmark,cardiac,stemness,celltype}.csv` — GSEA (Term, NES, NOM p-val, FDR q-val).
- `colocalization.csv` — every gene set × comparison (pooled, 4 cross- and 2 within-state replicate pairs): contact counts,
  per-state fold, hypergeometric p, Fisher odds ratio / p.
- `panel_table.csv` (5b), `supp_cardiac_sets_table.csv`, `supp_stemness_sets_table.csv` (5Sa c, d),
  `supp_coloc_vs_expression_table.csv` (5c, 5Sb) — pooled colocalization with the joint FDR, joined with GSEA.
- `supp_table_celltypes_perstate_and_fisher.csv` (5d), `supp_table_pathways_perstate_and_fisher.csv` (5c) — reference
  tables (per-state FDR + Fisher differential test), not plotted.
- `kcores/{sarcomere,cardiomyocyte}_{genes,states}.csv` + `kcore_summary.csv` (5e, 5f). Same-bin (cis) gene pairs have
  code 0 ("no interchromosomal contact").
- `replicate_agreement.csv` — Spearman agreement of each replicate pair with the pooled profile (reported in the
  figure-legend text of 5Sa b).

**FDR:** per-state enrichment is BH-corrected jointly over every pathway-level set tested (60-pathway panel (Hallmark +
custom) + 14 cardiac + 13 stemness = 87 sets), so a set has one q in every panel ("joint FDR"). Cell types (178) and the
Fisher tests are corrected separately.

**Re-reads:** the analysis builds some tables from CSVs read back from `data/` (and reads the GSEA tables with
`float_precision="round_trip"`), exactly as the published tables were built, so every value is bit-identical. pandas'
default CSV parser can differ from the in-memory float by up to ~2e-16 relative for folds / NES and up to ~9e-13 relative
for small p / q values — far below any reported precision. See the note in `psccm_analysis.ipynb` before changing a parser.

## Inputs
- Under `/data/land/3d-inter-paper-data/psc-cm-gse116862/`: LAS graphs `runs/pooled-v2/analysis/G_all_day00_vs_day80_filt_20.0.pkl`
  and `runs/replicate-resolved/analysis/G_all_<a>_vs_<b>_filt_20.0.pkl`; recount3 counts `rnaseq/recount3/gene_sums.G026.gz`.
- Gene sets: the 60-pathway panel as `../../metadata/hallmark_and_custom_pathways_unfiltered_aug_2025.json`
  (colocalization) and `.gmt` (GSEA; the same 60 sets); `../../metadata/psc-cm/` (main panel `psc_cm_pathways.json`,
  14 published cardiomyocyte sets `cardiomyocyte_published_sets.json`, 13 stemness sets `stemness_published_sets.json`,
  178 PanglaoDB cell types `panglaodb_celltypes.json`).
- `../../metadata/hg19/bins_df_protein_coding_genes_with_tss_hg19_250000.csv` (TSS genes and chromosome per graph node),
  `../../metadata/hg38/hg38_gtf.csv` (Ensembl ID → symbol).

## Notes
- Cleaned up 2026-10-21: `psccm_generate.ipynb` (colocalization) and the computations of `psccm.ipynb` (supplementary
  tables, k-cores, replicate agreement) moved into `psccm_analysis.ipynb`; the joint FDR / table building moved from
  `psccm_expression.ipynb` (which now only runs DESeq2 + GSEA) into the analysis; `psccm.ipynb` became the figures-only
  `psccm_figures.ipynb`. Supplementary panel files renamed to their final letters (`figS5{a-d}` → `fig5Sa_{a-d}`,
  `figS6{a,b}` → `fig5Sb_{a,b}`); `celltype_gsea.csv` → `gsea_celltype.csv`. The k-cores use the shared `../kcore.py`
  (same cores). All tables are unchanged.
- Legends framed in black (house style) on 2026-10-21: the 5b "FDR", 5f contact-state, 5Sa c / d and 5Sb b "GSEA FDR"
  legends (previously matplotlib's default grey frame). No other panel changed.
