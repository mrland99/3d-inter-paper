"""hg19 250 kb bins -> protein-coding genes by TSS (GENCODE v19), matching the muscle / CD8 (mm10) and fibroblast
(metadata/hg38/build_hg38_tss_annotation.py) convention
(each protein-coding gene assigned to the single bin containing its TSS).

Writes a table aligned row-for-row with bins_df_annotated_protein_coding_genes_hg19_250000.csv (row index = LAS graph
node id), column `protein_coding_genes_with_tss` = comma-separated gene names."""
import gzip
import re
import pandas as pd

HG19 = "/data/land/3d-inter-paper-data/hg19_data"
GTF = f"{HG19}/gencode.v19.annotation.gtf.gz"
BINS = f"{HG19}/bins_df_annotated_protein_coding_genes_hg19_250000.csv"     # existing bin grid (all-overlap genes)
OUT = "bins_df_protein_coding_genes_with_tss_hg19_250000.csv"   # written next to this script (run from metadata/hg19/)
RES = 250_000

# protein-coding gene TSSs (one per gene name; first occurrence in genomic order)
tss = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "gene" or 'gene_type "protein_coding"' not in f[8]:
            continue
        name = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
        chrom = f[0].replace("chr", "")
        pos = int(f[3]) if f[6] == "+" else int(f[4])
        tss.setdefault(name, (chrom, pos))
print(f"{len(tss)} protein-coding genes with a TSS")

bins = pd.read_csv(BINS, index_col=0)
bins["chrom"] = bins["chrom"].astype(str)
row_of = {(c, s): i for i, (c, s) in enumerate(zip(bins.chrom, bins.start))}
genes_in_bin = {i: [] for i in range(len(bins))}
unplaced = 0
for name, (chrom, pos) in tss.items():
    i = row_of.get((chrom, (pos - 1) // RES * RES))
    if i is None:
        unplaced += 1                                   # chrX/Y/M or contigs not in the autosomal grid
        continue
    genes_in_bin[i].append(name)

out = bins[["chrom", "start", "end"]].copy()
out["protein_coding_genes_with_tss"] = [",".join(genes_in_bin[i]) for i in range(len(bins))]
out.to_csv(OUT)
placed = sum(len(v) for v in genes_in_bin.values())
print(f"placed {placed} genes in {sum(1 for v in genes_in_bin.values() if v)} bins ({unplaced} not on the grid) -> {OUT}")
