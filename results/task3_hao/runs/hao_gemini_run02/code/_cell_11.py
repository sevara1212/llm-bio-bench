# Let's check marker genes for standard PBMC cell types!
markers = [
    'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', # T cells
    'IL7R', 'CCR7', 'S100A4', # CD4 Naive / Memory
    'NKG7', 'GNLY', 'NCAM1', 'FCGR3A', # NK / Cytotoxic
    'MS4A1', 'CD19', 'CD79A', 'CD79B', # B cells
    'CD14', 'LYZ', 'S100A9', 'S100A8', # CD14+ Monocytes
    'FCGR3A', 'MS4A7', # FCGR3A+ Monocytes
    'HLA-DRA', 'HLA-DQA1', 'CST3', # Antigen presenting / Monocytes / Dendritic
    'FCER1A', 'CST3', 'CLEC10A', # cDC
    'IL3RA', 'CLEC4C', 'LILRA4', 'TCF4', # pDC
    'PPBP', 'PF4', 'TUBB1', 'ITGA2B', # Platelets / Megakaryocytes
    'MZB1', 'SDC1', 'JCHAIN', # Plasma cells
    'HBA1', 'HBB' # Erythrocytes
]

import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, mask_var="highly_variable", svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42, key_added='leiden')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_genes = pd.DataFrame({group: [result['names'][group][i] for i in range(10)] for group in groups})
print("Top 10 marker genes per cluster (res=0.5):")
print(top_genes)