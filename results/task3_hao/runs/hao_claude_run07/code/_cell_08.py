import scanpy as sc
adata_hvg = sc.read_h5ad('adata_hvg.h5ad')

for res in [0.1,0.2,0.3,0.4]:
    sc.tl.leiden(adata_hvg, resolution=res, key_added=f'leiden_{res}')
    print(res, adata_hvg.obs[f'leiden_{res}'].nunique())

adata_hvg.write('adata_hvg.h5ad')