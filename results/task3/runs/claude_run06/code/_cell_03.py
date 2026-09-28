import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_mtx('filtered_gene_bc_matrices/hg19/matrix.mtx').T
genes = pd.read_csv('filtered_gene_bc_matrices/hg19/genes.tsv', sep='\t', header=None, names=['gene_id','gene_symbol'])
barcodes = pd.read_csv('filtered_gene_bc_matrices/hg19/barcodes.tsv', header=None, names=['barcode'])

adata.var_names = genes['gene_symbol'].values
adata.var['gene_id'] = genes['gene_id'].values
adata.obs_names = barcodes['barcode'].values
adata.var_names_make_unique()

print(adata)
adata.write('raw.h5ad')