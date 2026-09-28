with open('filtered_gene_bc_matrices/hg19/barcodes.tsv') as f:
    raw_bcs = [line.strip() for line in f]
print(f"Total raw barcodes: {len(raw_bcs)}")
print("Sample barcodes:", raw_bcs[:5])