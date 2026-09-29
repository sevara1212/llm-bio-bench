# Let's save the processed adata to disk so we don't need to recompute everything every time
adata.write_h5ad('processed.h5ad')
print("Saved processed.h5ad")