# Standard scanpy PBMC workflow:
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Filter cells with >2500 genes and >5% mitochondrial counts
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

print("Remaining cells after QC:", adata.n_obs)
print("Remaining genes:", adata.n_vars)

# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", adata.var.highly_variable.sum())

# Store raw
adata.raw = adata

# Regress out and scale
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.regress_out(adata_hvg, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata_hvg, max_value=10)

# PCA
sc.tl.pca(adata_hvg, svd_solver='arpack')
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

# Neighbors & UMAP & Leiden
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
print(adata.obs['leiden'].value_counts())