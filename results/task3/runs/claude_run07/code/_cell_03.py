import scanpy as sc
import pandas as pd
import numpy as np

sc.settings.verbosity = 1

adata = sc.read_mtx('filtered_gene_bc_matrices/hg19/matrix.mtx').T
genes = pd.read_csv('filtered_gene_bc_matrices/hg19/genes.tsv', header=None, sep='\t')
barcodes = pd.read_csv('filtered_gene_bc_matrices/hg19/barcodes.tsv', header=None, sep='\t')

adata.var['gene_ids'] = genes[0].values
adata.var_names = genes[1].values
adata.var_names_make_unique()
adata.obs_names = barcodes[0].values

print(adata)
adata.write('raw.h5ad')