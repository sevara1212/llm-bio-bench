import scanpy as sc
import pandas as pd
import numpy as np

# Load data
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
raw_barcodes = adata.obs_names.copy()

# QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter
filtered_idx = (adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5)
adata_filtered = adata[filtered_idx, :].copy()
sc.pp.filter_genes(adata_filtered, min_cells=3)

# Normalization & log
sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)

# HVG
sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_filtered.raw = adata_filtered

# Scale & PCA
adata_proc = adata_filtered[:, adata_filtered.var.highly_variable].copy()
sc.pp.scale(adata_proc, max_value=10)
sc.tl.pca(adata_proc, svd_solver='arpack')
sc.pp.neighbors(adata_proc, n_neighbors=10, n_pcs=40)

# Clustering at resolution 0.6
sc.tl.leiden(adata_proc, resolution=0.6, random_state=0)
adata_filtered.obs['leiden'] = adata_proc.obs['leiden']

# Let's inspect cluster markers
sc.tl.rank_genes_groups(adata_filtered, 'leiden', method='wilcoxon')
result = adata_filtered.uns['rank_genes_groups']

# Let's evaluate key genes
genes_to_check = ['CD3D', 'IL7R', 'CD8A', 'NKG7', 'GNLY', 'CD14', 'LYZ', 'FCGR3A', 'MS4A1', 'FCER1A', 'CST3', 'PPBP']
cluster_profiles = {}
for cl in sorted(adata_filtered.obs['leiden'].unique()):
    top_markers = [result['names'][cl][i] for i in range(5)]
    exprs = {}
    for g in genes_to_check:
        if g in adata_filtered.raw.var_names:
            exprs[g] = float(np.asarray(adata_filtered.raw[adata_filtered.obs['leiden'] == cl, g].X.todense()).mean())
    cluster_profiles[cl] = (top_markers, exprs)
    print(f"Cluster {cl} top markers: {top_markers}")
    print(f"  Exprs: {', '.join([f'{k}: {v:.2f}' for k, v in exprs.items() if v > 0.3])}")

# Mapping:
# 0: CD4 T cells / T cells (CD3D high, IL7R high, ribosomal, LDHB)
# 1: CD14+ Monocytes (LYZ high, S100A8/9 high, CD14 high)
# 2: NK cells (NKG7 high, GNLY high, CST7 high) or CD8 T cells / NK cells
# 3: B cells (CD79A/B high, MS4A1 high, HLA high)
# 4: FCGR3A+ Monocytes (FCGR3A high, LST1 high)
# 5: Dendritic cells (FCER1A high, CST3 high, HLA high)
# 6: Megakaryocytes (PPBP high, PF4 high)

cluster_labels = {
    '0': 'T cells',
    '1': 'CD14+ Monocytes',
    '2': 'NK cells',
    '3': 'B cells',
    '4': 'FCGR3A+ Monocytes',
    '5': 'Dendritic cells',
    '6': 'Megakaryocytes'
}

# Assign to filtered cells
adata_filtered.obs['cell_type'] = adata_filtered.obs['leiden'].map(cluster_labels)

# For any cells filtered out (if full barcodes needed), let's see:
df = pd.DataFrame({'barcode': adata_filtered.obs_names, 'cell_type': adata_filtered.obs['cell_type'].values})
df.to_csv('labels.csv', index=False)
print("\nSaved labels.csv with shape:", df.shape)
print(df['cell_type'].value_counts())
print("\nFirst 5 rows:")
print(df.head())