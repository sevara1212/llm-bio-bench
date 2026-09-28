import scanpy as sc

# The canonical pbmc3k tutorial sets:
# sc.tl.pca(adata, svd_solver='arpack')
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
# sc.tl.umap(adata)
# sc.tl.leiden(adata) # default resolution is 1.0 or 0.8
adata = sc.read_h5ad('adata_processed.h5ad')
for res in [0.7, 0.8, 0.9, 1.0, 1.2]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    print(f"Res {res}: {len(adata.obs[f'leiden_{res}'].unique())} clusters")