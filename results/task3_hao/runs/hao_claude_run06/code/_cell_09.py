import scanpy as sc

adata = sc.read_h5ad('clustered.h5ad')

# Try a couple of resolutions to find a reasonable number of clusters for cell typing
for res in [0.3, 0.4, 0.5, 0.6, 0.8]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    print(res, adata.obs[f'leiden_{res}'].nunique())