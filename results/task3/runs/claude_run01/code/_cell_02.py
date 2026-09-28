import scanpy as sc
adata = sc.read_mtx('filtered_gene_bc_matrices/hg19/matrix.mtx').T
print(adata.shape)
import pandas as pd
genes = pd.read_csv('filtered_gene_bc_matrices/hg19/genes.tsv', sep='\t', header=None)
barcodes = pd.read_csv('filtered_gene_bc_matrices/hg19/barcodes.tsv', sep='\t', header=None)
print(genes.shape, barcodes.shape)
print(genes.head())
print(barcodes.head())