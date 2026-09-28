import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
raw_barcodes = pd.Series(index=adata.obs_names, dtype=str)

# QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Standard filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

# Normalization and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# HVGs
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_sub = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_sub, max_value=10)
sc.tl.pca(adata_sub, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_sub, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.umap(adata_sub, random_state=0)

# Check resolution 0.5 clusters
sc.tl.leiden(adata_sub, resolution=0.5, random_state=0)
adata.obs['cluster'] = adata_sub.obs['leiden']

# Let's inspect marker expression across the 6 clusters:
cluster_names = {
    '0': 'T cells',
    '1': 'Monocytes',
    '2': 'NK cells',
    '3': 'B cells',
    '4': 'Dendritic cells',
    '5': 'Platelets'
}
adata.obs['cell_type'] = adata.obs['cluster'].map(cluster_names)

# Check markers
for gene in ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD14', 'FCGR3A', 'MS4A7', 'NKG7', 'GNLY', 'MS4A1', 'CD79A', 'FCER1A', 'CST3', 'PPBP', 'PF4']:
    if gene in adata.raw.var_names:
        print(f"\n{gene}:")
        for cl in sorted(adata.obs['cluster'].unique(), key=lambda x: int(x)):
            val = adata.raw[adata.obs['cluster'] == cl, gene].X.mean()
            print(f"  Cl {cl} ({cluster_names[cl]}): {val:.3f}")