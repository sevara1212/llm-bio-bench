import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
raw_barcodes = pd.read_csv('filtered_gene_bc_matrices/hg19/barcodes.tsv', header=None, sep='\t')[0].values
print("Total raw barcodes:", len(raw_barcodes))