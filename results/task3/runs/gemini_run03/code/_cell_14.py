import scanpy as sc
import numpy as np
import pandas as pd

# 1. Load data
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
print(f"Initial cells: {adata.n_obs}, genes: {adata.n_vars}")

# 2. QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Standard filtering criteria for PBMC 3k:
# Filter cells with < 200 genes, > 2500 genes, or > 5% mitochondrial counts
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[(adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5), :].copy()
print(f"Filtered cells: {adata.n_obs}, genes: {adata.n_vars}")

# 3. Normalization and log transformation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# 4. Feature selection and scaling
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_sub = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_sub, max_value=10)

# 5. Dimensionality reduction and Clustering
sc.tl.pca(adata_sub, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_sub, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.umap(adata_sub, random_state=0)
sc.tl.leiden(adata_sub, resolution=0.5, random_state=0)

adata.obs['leiden'] = adata_sub.obs['leiden']

# 6. Marker gene identification
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# 7. Cell type annotation
cluster_annotations = {
    '0': 'T cells',
    '1': 'Monocytes',
    '2': 'NK cells',
    '3': 'B cells',
    '4': 'Dendritic cells',
    '5': 'Platelets'
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)

# 8. Save labels.csv
df = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type']
})
df.to_csv('labels.csv', index=False)
print("Saved labels.csv with shape:", df.shape)
print("Cell type distribution:")
print(df['cell_type'].value_counts())