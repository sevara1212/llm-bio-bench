import scanpy as sc

adata = sc.read_h5ad('raw_counts.h5ad')
print("Barcodes head:")
print(adata.obs_names[:10].tolist())

# Check gene counts across cells
sc.pp.filter_genes(adata, min_cells=3)
print("Genes remaining after min_cells=3:", adata.n_vars)