# Figure 4 — interchromosomal colocalization in aging muscle stem cells

Interchromosomal (*trans*) LAS pathway colocalization in mouse skeletal-muscle stem cells at 250 kb: Zhao et al. 2023
(main figure; young 2 mo / aged 23–24 mo / geriatric 28–32 mo) and the independent NIA replication (Shi et al. 2025;
Hi-C GSE293131, RNA-seq GSE293130; young 2–4 / old 20–24 / geriatric 27–32 mo, male).

## Notebooks (run in this order, from this folder)
Zhao (Figure 4, Supplementary 4Sa, 4Sb):
1. `muscle_analysis.ipynb` — every computation → `data/muscle/`: per-state hypergeometric colocalization, Fisher
   differential specificity and fold enrichment for the three contrasts; GSEA (shared `../gsea.py`); the locus-matched
   permutation test (shared `../locus_null.py`) with its per-age expression strata; the geriatric-specific k-cores (shared
   `../kcore.py`); the age trajectory (excess over a random-gene null, net-accumulation sweep). ~15 min (the
   net-accumulation sweep ~10 min).
2. `muscle_figures.ipynb` — reads only `data/muscle/` and draws every panel → `figures/<panel>.svg`, one `###` section
   per final panel. ~30 s.
3. `muscle_pages.ipynb` — assembles `Figure4Sb`, cropped tight → SVG/PDF/PNG in `figures/`, page SVG in
   `../../illustrator-figures/`. Figure 4 and 4Sa are assembled in Illustrator (`../../illustrator-figures/Figure4.ai`,
   `Figure4Sa.ai`; panel 4a, the schematic, lives there), with the panel SVGs embedded.

NIA replication (Supplementary 4S NIA), independent of the above:
4. `nia_replication.ipynb` — self-contained: analysis → `data/nia/`, panels a–f (`###` sections for a, b / c, d /
   e, f, one per contrast pair) → `figures/fig4S_NIA_*.svg`, and the page `Figure4S_NIA` (SVG/PDF/PNG in `figures/`, SVG
   in `../../illustrator-figures/`). ~3 min.

To run all notebooks headless, e.g. with `nbclient`:

```python
import nbformat
from nbclient import NotebookClient
KERNEL = "python3"                                       # any Jupyter kernel with the required packages
for nb_name in ["muscle_analysis.ipynb", "muscle_figures.ipynb", "muscle_pages.ipynb", "nia_replication.ipynb"]:
    nb = nbformat.read(nb_name, as_version=4)
    NotebookClient(nb, kernel_name=KERNEL, timeout=3600,
                   resources={"metadata": {"path": "."}}).execute()   # run from this folder so data/ and figures/ resolve
    nbformat.write(nb, nb_name)
```

## Panel → notebook section → file
Each Zhao panel has a `###` section in `muscle_figures.ipynb` under the `##` section below; its computation is under the
same-named `##` section of `muscle_analysis.ipynb` (4c combines the colocalization and GSEA sections). The NIA panels are
`###` sections in `nia_replication.ipynb`, computed in its sections of the same names.

| Panel | Content | `##` section | Panel file (`figures/`) |
|---|---|---|---|
| Figure 4a | experiment schematic | — (Illustrator) | — |
| Figure 4b | colocalized pathways, young vs geriatric | Pathway colocalization | `fig4b_colocalization_young_vs_geriatric` |
| Figure 4c | colocalization vs expression, young vs geriatric | Colocalization vs expression (from Pathway colocalization + Pathway expression) | `fig4c_concordance_young_vs_geriatric` |
| Figure 4d | IFN-γ response geriatric-specific core | State-specific k-cores | `fig4d_ifng_core` |
| Figure 4e | colocalized pathways, aged vs geriatric | Pathway colocalization | `fig4e_colocalization_aged_vs_geriatric` |
| Figure 4f | colocalization over age, relative to young | Age trajectory of colocalization | `fig4f_colocalization_over_age` |
| Figure 4g | net age-gained colocalization vs random genes | Age trajectory of colocalization | `fig4g_net_accumulation` |
| Supplementary 4Sa, panel a | GSEA, young vs geriatric | Pathway expression (GSEA) | `fig4Sa_a_gsea_young_vs_geriatric` |
| Supplementary 4Sa, panel b | colocalized pathways, young vs aged | Pathway colocalization | `fig4Sa_b_colocalization_young_vs_aged` |
| Supplementary 4Sa, panel c | GSEA, young vs aged | Pathway expression (GSEA) | `fig4Sa_c_gsea_young_vs_aged` |
| Supplementary 4Sa, panel d | GSEA, aged vs geriatric | Pathway expression (GSEA) | `fig4Sa_d_gsea_aged_vs_geriatric` |
| Supplementary 4Sa, panel e | absolute colocalization level over age | Age trajectory of colocalization | `fig4Sa_e_absolute_colocalization` |
| Supplementary 4Sb, panels a, b | locus-matched test, young vs geriatric / aged vs geriatric | Locus-matched permutation test | `fig4Sb_{a,b}_locus_matched_*` |
| Supplementary 4Sb, panel c | IFN-α and complement geriatric-specific cores | State-specific k-cores | `fig4Sb_c_ifna_complement_cores` |
| Supplementary 4S NIA, panels a, b | NIA colocalized pathways (young vs geriatric / old) | Pathway colocalization (`nia_replication.ipynb`) | `fig4S_NIA_{a,b}_colocalization_*` |
| Supplementary 4S NIA, panels c, d | NIA locus-matched test | Locus-matched permutation test (`nia_replication.ipynb`) | `fig4S_NIA_{c,d}_locus_matched_*` |
| Supplementary 4S NIA, panels e, f | NIA colocalization vs expression | Pathway colocalization + Pathway expression (`nia_replication.ipynb`) | `fig4S_NIA_{e,f}_concordance_*` |

## Outputs in `data/`
- `muscle/` (Zhao): `pathway_enrichment_<pair>_filt_20.0.csv` (per-state hypergeometric, BH over the 60 sets),
  `differential_specificity_<pair>_filt_20.0.csv` (Fisher test of the state skew; reference only, not plotted),
  `fold_enrichment_<pair>.csv`, `gsea_<pair>.csv`, `expression_by_age.csv` (per-age mean salmon counts, the expression
  strata of the locus-matched test), `locus_matched_{young,aged}_vs_geriatric.csv`,
  `kcores/{ifng,ifna,complement}_{genes,states}.csv` + `kcore_summary.csv`, `trajectory_{relative,absolute}.csv`,
  `net_accumulation_sweep.csv`.
- `nia/`: `pathway_enrichment_young_vs_<older>.csv` (incl. fold), `gsea_young_vs_<older>.csv`,
  `locus_matched_young_vs_<older>.csv`.

## Inputs
- Under `/data/land/3d-inter-paper-data/`: Zhao `muscle-sc-zhao/runs/pooled-v2/analysis/G_all_<pair>_filt_20.0.pkl`
  (age-specific contact graphs), per-age LAS coolers
  `muscle-sc-zhao/runs/pooled-v2/cool/<age>_INTER_LAS_250000_filt_20.0.cool`, salmon quantifications `muscle-sc-zhao/raw-data/rna/quant/<SRR>/quant.sf` and PyDESeq2 results
  `muscle-sc-zhao/raw-data/rna/de_<pair>.csv`; NIA graphs `muscle-sc-nia/runs/pooled-v2/analysis/`
  and PyDESeq2 results `muscle-sc-nia/rna/de_young_vs_<older>.csv` (male samples, + = up in the older group); mm10 250 kb
  bins `cd8-hu/annot/mm10_250kb/mm10_250kb_bins_annotated.csv` (row = graph node / cooler bin id) and the transcript →
  gene table `cd8-hu/annot/gencode.vM25.tx2gene.tsv`.
- In `data/nia/`: `expression_by_age.csv`, the per-age mean DESeq2-normalized counts of the 12 male NIA samples (size
  factors fit on all 12), made with the NIA PyDESeq2 results by the same RNA-seq processing step; it sets the expression
  strata of the NIA locus-matched test.
- Shared: `../../metadata/mm10_hallmark_and_custom_pathways.json` (60 Hallmark + custom sets, mouse).

## Notes
- Cleaned up 2026-10-20: `muscle_generate.ipynb`, `muscle.ipynb`, `age_progressive_colocalization.ipynb`,
  `figure4Sb.ipynb` and `nia_supplement.ipynb` were merged into analysis / figures / pages notebooks; panel files were
  renamed to their final letters; the exploratory trajectory plots (raw counts, fold, null band, validation, sweep volcano)
  and the unused aged and aged-vs-geriatric concordance scatters were dropped; `all_pathways_sweep.csv` moved to
  `data/muscle/net_accumulation_sweep.csv`. The k-cores now use the shared `../kcore.py` (same cores).
- The NIA locus-matched test (4S NIA c, d) is now run here with the shared `../locus_null.py`, as for Zhao (4Sb a, b);
  it was previously copied from `workspace/coloc-null/coloc_null_general.py` (variant `matched_max`).
- The published panels were made with matplotlib 3.11; older versions (e.g. 3.5) render the SVGs slightly differently.
- The NIA locus-matched p-values differ from the previous copy by Monte Carlo noise only (same method, inputs and seeds;
  the shared module visits strata and loci in a fixed order); one set crosses p = 0.05 (IL6 JAK STAT3 signaling, young,
  young vs old: 0.049 → 0.053).
