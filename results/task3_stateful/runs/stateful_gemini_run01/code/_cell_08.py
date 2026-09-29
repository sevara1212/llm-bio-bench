# Let's perform standard scanpy PBMC3k workflow!
# 1. Quality Control
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Standard filtering in scanpy tutorial:
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

print("Remaining cells:", adata.n_obs, "genes:", adata.n_vars)

# 2. Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save raw counts / state
adata.raw = adata

# 3. Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", adata.var.highly_variable.sum())

# 4. Scale and PCA
adata_processed = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_processed, max_value=10)
sc.tl.pca(adata_processed, svd_solver='arpack', random_state=0)

# Add PCA back or compute neighbors on adata_processed
sc.pp.neighbors(adata_processed, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.leiden(adata_processed, resolution=0.5, random_state=0)
adata.obs['leiden'] = adata_processed.obs['leiden']

print("Cluster counts:")
print(adata.obs['leiden'].value_counts())