# Let's check markers across different resolutions or higher resolution (like 0.8 or 1.0)
for res in [0.4, 0.5, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata_hvg, resolution=res, random_state=0, key_added=f'leiden_{res}')
    print(f"Res {res}: {adata_hvg.obs[f'leiden_{res}'].nunique()} clusters: {adata_hvg.obs[f'leiden_{res}'].value_counts().to_dict()}")