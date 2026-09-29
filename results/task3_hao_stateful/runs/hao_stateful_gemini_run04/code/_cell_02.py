adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("obs head:", adata.obs.head())
print("var head:", adata.var.head())
print("X shape:", adata.X.shape)
if hasattr(adata.X, 'toarray'):
    print("X min/max:", adata.X.min(), adata.X.max())
else:
    print("X min/max:", np.min(adata.X), np.max(adata.X))