# Let's perform standard Scanpy PBMC workflow:
# 1. Filter genes: min_cells=3
# 2. Total counts normalisation to 10,000 (target_sum=1e4)
# 3. log1p
# 4. Highly variable genes (e.g. flavor='seurat' or 'cell_ranger', n_top_genes=2000)
# 5. Scale data (and save raw/log1p for marker genes)
# 6. PCA (e.g. n_comps=30 or 50)
# 7. Neighbors (n_neighbors=10 or 15, n_pcs=20 or 30)
# 8. Leiden / Louvain clustering
# 9. Find marker genes: rank_genes_groups
# 10. Check markers for each cluster and assign cell types

# Let's test this in Python
adata_proc = adata.copy()
sc.pp.filter_genes(adata_proc, min_cells=3)
sc.pp.normalize_total(adata_proc, target_sum=1e4)
sc.pp.log1p(adata_proc)
adata_proc.raw = adata_proc  # save raw for marker gene test

sc.pp.highly_variable_genes(adata_proc, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(f"Highly variable genes: {adata_proc.var.highly_variable.sum()}")

adata_proc_hvg = adata_proc[:, adata_proc.var.highly_variable].copy()
sc.pp.scale(adata_proc_hvg, max_value=10)
sc.tl.pca(adata_proc_hvg, svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata_proc_hvg, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata_proc_hvg, resolution=0.5, random_state=42)

adata_proc.obs['leiden_0.5'] = adata_proc_hvg.obs['leiden']
sc.tl.leiden(adata_proc_hvg, resolution=0.8, random_state=42)
adata_proc.obs['leiden_0.8'] = adata_proc_hvg.obs['leiden']

print("Clusters at res 0.5:", adata_proc.obs['leiden_0.5'].value_counts())
print("Clusters at res 0.8:", adata_proc.obs['leiden_0.8'].value_counts())