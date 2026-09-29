import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=3)

# Keep raw counts in adata.raw
adata.raw = adata

# Normalisation and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Scale & PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata, resolution=0.3, key_added='leiden', random_state=42)

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
markers = adata.uns['rank_genes_groups']

canonical = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'IL7R', 'GNLY', 'NKG7', 'MS4A1', 'CD14', 'FCGR3A', 'FCER1A', 'PPBP', 'LILRA4']

for g in sorted(adata.obs['leiden'].unique(), key=int):
    top_m = [markers['names'][g][i] for i in range(8)]
    # mean expr of canonical markers in raw
    cells = adata.obs['leiden'] == g
    mean_canon = {c: f"{adata.raw[cells, c].X.mean():.2f}" for c in canonical if c in adata.raw.var_names}
    print(f"Cluster {g} (n={cells.sum()}): {', '.join(top_m)}")
    print(f"   Key markers: CD3D={mean_canon.get('CD3D')}, CD4={mean_canon.get('CD4')}, CD8A={mean_canon.get('CD8A')}, GNLY={mean_canon.get('GNLY')}, NKG7={mean_canon.get('NKG7')}, MS4A1={mean_canon.get('MS4A1')}, CD14={mean_canon.get('CD14')}, FCGR3A={mean_canon.get('FCGR3A')}, FCER1A={mean_canon.get('FCER1A')}, PPBP={mean_canon.get('PPBP')}, LILRA4={mean_canon.get('LILRA4')}")