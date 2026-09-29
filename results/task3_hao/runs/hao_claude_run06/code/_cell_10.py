import scanpy as sc

adata = sc.read_h5ad('clustered.h5ad')

for res in [0.1, 0.15, 0.2, 0.25]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    print(res, adata.obs[f'leiden_{res}'].nunique(), adata.obs[f'leiden_{res}'].value_counts().to_dict())