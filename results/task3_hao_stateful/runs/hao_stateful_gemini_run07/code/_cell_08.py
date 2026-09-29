# Check if there are uns, layers, or other metadata in raw_counts.h5ad
print("adata.uns:", adata.uns.keys())
print("adata.layers:", adata.layers.keys())
print("adata.obs:", adata.obs.columns.tolist())