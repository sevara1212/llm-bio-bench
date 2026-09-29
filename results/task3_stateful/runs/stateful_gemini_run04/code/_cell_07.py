# Let's run the tutorial steps and check clusters
# Copy adata
adata_filtered = adata.copy()
sc.pp.filter_cells(adata_filtered, min_genes=200)
sc.pp.filter_genes(adata_filtered, min_cells=3)
adata_filtered = adata_filtered[adata_filtered.obs.n_genes_by_counts < 2500, :]
adata_filtered = adata_filtered[adata_filtered.obs.pct_counts_mt < 5, :].copy()

# Normalization & Log
sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)

# Highly variable genes
sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_filtered.raw = adata_filtered

adata_proc = adata_filtered[:, adata_filtered.var.highly_variable].copy()
sc.pp.scale(adata_proc, max_value=10)
sc.tl.pca(adata_proc, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_proc, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.leiden(adata_proc, resolution=0.5, random_state=0)

adata_filtered.obs['leiden'] = adata_proc.obs['leiden']
print(adata_filtered.obs['leiden'].value_counts())