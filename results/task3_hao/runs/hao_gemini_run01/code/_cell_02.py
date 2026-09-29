import scanpy as sc
import anndata as ad

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("obs head:", adata.obs.head())
print("var head:", adata.var.head())
print("X shape:", adata.X.shape)
print("X type:", type(adata.X))
print("min/max X:", adata.X.min(), adata.X.max())