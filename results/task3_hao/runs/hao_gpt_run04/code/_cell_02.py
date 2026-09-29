import scanpy as sc
adata=sc.read_h5ad('raw_counts.h5ad')
print(adata)
print(adata.obs.head(),adata.var.head(),sep='\n')
print(adata.X.min(),adata.X.max(), adata.X.sum())
print(adata.var_names[:10].tolist())
print(adata.obs.columns.tolist(),adata.var.columns.tolist())