for res in [0.4, 0.5, 0.6, 0.7, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}', random_state=0)
    print(f"Res {res}: {adata.obs[f'leiden_{res}'].nunique()} clusters")
    print(adata.obs[f'leiden_{res}'].value_counts().to_dict())