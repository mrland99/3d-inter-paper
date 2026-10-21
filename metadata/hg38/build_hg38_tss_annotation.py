"""hg38 250 kb bins -> protein-coding genes by TSS, one bin per gene (same rule as metadata/hg19/build_hg19_tss_annotation.py
and the mm10 annotation used by Figs 4 and 6).

TSS = 5' end of the GENCODE `gene` record: start on the + strand, end on the - strand (the gene record spans all of its
transcripts, so this is the most upstream transcript start). Row index = hg38 250 kb bin = fibroblast LAS graph node id.
The fibroblast graph's own gene lists (`protein_coding_genes_with_tss`, Jana's annotation) hold every bin a gene body
overlaps; the notebooks replace them with this table at load. Writes bins_protein_coding_genes_tss_hg38_250000.csv;
run from metadata/hg38/."""
import pandas as pd

RES = 250_000

gtf = pd.read_csv("hg38_gtf.csv", index_col=0, usecols=["Unnamed: 0", "chr", "feature", "start", "end", "strand",
                                                         "gene_name", "gene_type"])
genes = gtf[(gtf.feature == "gene") & (gtf.gene_type == "protein_coding")]
genes = genes.drop_duplicates("gene_name")                       # one TSS per gene name: first record in genomic order
tss = genes.start.where(genes.strand == "+", genes.end)          # 1-based TSS
print(f"{len(genes)} protein-coding genes with a TSS")

bins = pd.read_csv("250kb_loci_annotations/bins_df_annotated_protein_coding_genes_and_DEGs.csv", index_col=0)
row_of = {(str(c), s): i for i, (c, s) in enumerate(zip(bins.chrom, bins.start))}
genes_in_bin = {i: [] for i in range(len(bins))}
unplaced = 0
for name, chrom, pos in zip(genes.gene_name, genes.chr, tss):
    i = row_of.get((chrom.removeprefix("chr"), (int(pos) - 1) // RES * RES))
    if i is None:
        unplaced += 1                                            # chrX/Y/M and contigs are not on the autosomal grid
        continue
    genes_in_bin[i].append(name)
print(f"placed {len(genes) - unplaced}, unplaced {unplaced}")

out = bins[["chrom", "start", "end"]].copy()
out["protein_coding_genes_tss"] = [",".join(genes_in_bin[i]) for i in range(len(bins))]
out.to_csv("bins_protein_coding_genes_tss_hg38_250000.csv")
print("wrote bins_protein_coding_genes_tss_hg38_250000.csv")
