import scanpy as sc
import anndata as ad

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("X shape:", adata.X.shape)
print("obs head:\n", adata.obs.head())
print("var head:\n", adata.var.head())