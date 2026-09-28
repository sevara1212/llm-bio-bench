# Let's check how many barcodes are in barcodes.tsv
with open('filtered_gene_bc_matrices/hg19/barcodes.tsv') as f:
    raw_bcs = [l.strip() for l in f]
print("Total barcodes in raw matrix:", len(raw_bcs))