import scanpy as sc
adata = sc.read_mtx('filtered_gene_bc_matrices/hg19/matrix.mtx').T
genes = [l.strip().split('\t') for l in open('filtered_gene_bc_matrices/hg19/genes.tsv')]
barcodes = [l.strip() for l in open('filtered_gene_bc_matrices/hg19/barcodes.tsv')]
print(adata.shape, len(genes), len(barcodes))
print(genes[:5])
print(barcodes[:5])