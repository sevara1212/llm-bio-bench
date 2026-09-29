# Check how many barcodes are in the raw barcodes.tsv
raw_barcodes = pd.read_csv('filtered_gene_bc_matrices/hg19/barcodes.tsv', header=None)[0].values
print(f"Total raw barcodes: {len(raw_barcodes)}")