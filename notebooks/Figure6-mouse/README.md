# Figure 6 — interchromosomal colocalization of CD8 T-cell exhaustion programs (Hu + Li)

Interchromosomal (*trans*) LAS pathway colocalization at 250 kb in progenitor-exhausted (Tpex) vs effector-exhausted
(Tex) CD8 T cells, in two independent datasets: Hu et al. 2026 (LCMV clone 13 chronic infection) and Li et al. 2025
(E.G7-OVA tumor).

## Notebooks (run in this order, from this folder)
1. `cd8_analysis.ipynb` — every computation → `data/`: per-state hypergeometric colocalization and the Fisher
   differential test (both datasets), contact composition of the five program sets, GSEA of the 17 CD8-specific sets on
   the Hu ranking (shared `../gsea.py`), the all-17 CD8-specific set table, the degree-preserving null, the gene overlap
   of the five sets, and the Tex-specific effector k-cores with their Hu/Li overlap test (shared `../kcore.py`). ~85 s.
2. `cd8_figures.ipynb` — reads only `data/` and draws the panels → `figures/<panel>.svg`, one `###` section per final
   panel. ~5 s.

Figure 6 is assembled in Illustrator (`../../illustrator-figures/Figure6.ai`; panel 6a, the schematic, lives there), with
the panel SVGs embedded, so there is no pages notebook.

To run both notebooks headless, e.g. with `nbclient`:

```python
import nbformat
from nbclient import NotebookClient
KERNEL = "python3"                                       # any Jupyter kernel with the required packages
for nb_name in ["cd8_analysis.ipynb", "cd8_figures.ipynb"]:
    nb = nbformat.read(nb_name, as_version=4)
    NotebookClient(nb, kernel_name=KERNEL, timeout=1800,
                   resources={"metadata": {"path": "."}}).execute()   # run from this folder so data/ and figures/ resolve
    nbformat.write(nb, nb_name)
```

## Panel → notebook section → file
Each panel has a `###` section in `cd8_figures.ipynb`; its computation is under the listed `##` section(s) of
`cd8_analysis.ipynb`.

| Panel | Content | Analysis `##` section | Panel file (`figures/`) |
|---|---|---|---|
| Figure 6a | datasets and cell states (schematic) | — (Illustrator) | — |
| Figure 6b | GSEA of the five program sets, Tex vs Tpex, Hu (LCMV infection) | Pathway expression | `fig6b_gsea_hu` |
| Figure 6c | colocalization of the five program sets, Hu and Li | Pathway colocalization (set shares, FDR) + Contact composition (genome-wide shares) | `fig6c_colocalization_hu_li` |
| Figure 6d | Tex-specific effector k-cores, Hu and Li | Effector k-cores | `fig6d_effector_cores` |

## Supplementary table → notebook section → file (`data/`)

| Table | Analysis `##` section | File |
|---|---|---|
| Fisher differential test (Tex vs Tpex skew per set) | Pathway colocalization | `{cd8-hu,li-2025}/differential_specificity_filt_20.0.csv` |
| Contact composition of the five program sets | Contact composition | `{cd8-hu,li-2025}/contact_composition_5sets.csv` |
| GSEA of the 17 CD8-specific sets (Hu) | Pathway expression | `cd8-hu/gsea_hu_tex_vs_tpex.csv` |
| All 17 CD8-specific sets: fold, FDR (Hu, Li), GSEA | CD8-specific gene sets | `cd8_specific_sets_colocalization.csv` |
| Degree-preserving null | Degree-preserving null | `{cd8-hu,li-2025}/degree_null_5sets.csv` |
| Hu/Li effector-core overlap | Effector k-cores | `kcores/core_overlap.csv` |

## Outputs in `data/`
- Per dataset: `{cd8-hu,li-2025}/pathway_enrichment_filt_20.0.csv` (per-state hypergeometric; BH over the 77
  finalized gene sets), `differential_specificity_filt_20.0.csv` (Fisher; BH over the 77 finalized gene sets),
  `contact_composition_5sets.csv`, `degree_null_5sets.csv`.
- Also: `program_sets.csv` (the five program sets, ids / clean
  names / poles, read by the figures notebook), `cd8-hu/gsea_hu_tex_vs_tpex.csv` (FDR over the 17 CD8-specific sets),
  `cd8_specific_sets_colocalization.csv`, `kcores/{hu,li}_effector_{genes,states}.csv`, `kcores/kcore_summary.csv`,
  `kcores/core_overlap.csv`, `overlap_5sets.csv` (pairwise shared genes of the five sets; QC).

## Inputs (under `/data/land/3d-inter-paper-data/`)
- `cd8-hu/runs/pooled-v2/analysis/G_all_filt_20.0.pkl`, `li-2025/runs/pooled-v2/analysis/G_all_filt_20.0.pkl` —
  cell-state-specific contact graphs (edge types `T_PEX`, `T_EX_EFF`, `Shared`).
- `cd8-hu/annot/mm10_250kb/`: `mm10_250kb_bins_annotated.csv` (row = node id), `mm10_pathways_final.json` (77 finalized
  gene sets: 60 unbiased + 17 CD8-specific; the multiple-testing family),
  `pathway_catalog.csv` (clean name, tier, source, pole per finalized set), `hu_tex_vs_tpex_s2n.rnk` (Hu bulk RNA-seq
  Tex-vs-Tpex signal-to-noise ranking over the ~3,800-gene DEG list published as a supplementary table).

## Notes
- All colocalization FDRs (per-state hypergeometric, Fisher) are BH over the 77 finalized gene sets; the GSEA FDR is over
  the 17 CD8-specific sets.
- The panels were made with matplotlib 3.11.
