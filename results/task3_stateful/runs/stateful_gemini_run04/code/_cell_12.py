# Let's inspect resolutions for Leiden: e.g. 0.5, 0.6, 0.7, 0.8, 1.0
for res in [0.4, 0.5, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata_proc, resolution=res, random_state=0)
    print(f"Res {res}: n_clusters = {adata_proc.obs['leiden'].nunique()}")