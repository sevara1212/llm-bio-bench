# In scanpy PBMC tutorial:
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
# sc.tl.umap(adata)
# sc.tl.leiden(adata) # default resolution is 1.0 or 0.8
# Let's check with default resolution or resolution 0.8 / 0.5 / etc.
for res in [0.4, 0.5, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata_processed, resolution=res, random_state=0)
    print(f"res={res}: {len(adata_processed.obs['leiden'].unique())} clusters: {adata_processed.obs['leiden'].value_counts().to_dict()}")