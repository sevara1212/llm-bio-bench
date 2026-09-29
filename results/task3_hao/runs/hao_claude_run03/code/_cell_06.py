import scanpy as sc

adata = sc.read_h5ad('qc_filtered.h5ad')

# save raw counts layer
adata.layers['counts'] = adata.X.copy()

# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
print(adata.var['highly_variable'].sum())

adata.raw = adata

adata_hvg = adata[:, adata.var.highly_variable].copy()

sc.pp.scale(adata_hvg, max_value=10)

sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)

adata_hvg.write('pca.h5ad')