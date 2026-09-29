# Let's perform standard Scanpy PBMC workflow:
# 1. Filter genes: e.g. min_cells=3
# 2. Normalize to 10,000 (target_sum=1e4)
# 3. log1p
# 4. Highly variable genes (e.g. min_mean=0.0125, max_mean=3, min_disp=0.5 or n_top_genes=2000)
# 5. Scale data, regress out or PCA
# 6. Compute PCA (n_comps=30 or 50)
# 7. Compute neighborhood graph (sc.pp.neighbors)
# 8. Run UMAP and Leiden clustering (resolution=0.5 or 0.8)
# Let's test this in Python and check cluster markers.

adata_proc = adata.copy()
sc.pp.filter_genes(adata_proc, min_cells=3)
adata_proc.raw = adata_proc # store raw counts before normalization? Or after log1p?
# Usually sc.pp.normalize_total -> sc.pp.log1p -> adata_proc.raw = adata_proc
sc.pp.normalize_total(adata_proc, target_sum=1e4)
sc.pp.log1p(adata_proc)
adata_proc.raw = adata_proc

sc.pp.highly_variable_genes(adata_proc, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", adata_proc.var['highly_variable'].sum())

adata_hvg = adata_proc[:, adata_proc.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata_hvg, random_state=42)
sc.tl.leiden(adata_hvg, resolution=0.5, random_state=42)

adata_proc.obs['leiden_0.5'] = adata_hvg.obs['leiden']
print(adata_proc.obs['leiden_0.5'].value_counts())