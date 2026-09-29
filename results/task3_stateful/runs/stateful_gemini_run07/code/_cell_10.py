# Let's perform the standard Scanpy pipeline step by step:
# 1. Filter cells and genes
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Standard filtering in scanpy pbmc3k tutorial:
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

# Normalization & log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Save raw
adata.raw = adata

# Filter to HVG and scale
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver='arpack')

# Neighbors & UMAP & Leiden
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.leiden(adata, resolution=0.5, random_state=0)
print(adata.obs['leiden'].value_counts())