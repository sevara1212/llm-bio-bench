import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_filtered.h5ad')
genes = ['CD4', 'IL7R', 'CCR7', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'CD14', 'FCGR3A', 'MS4A1', 'FCER1A', 'PPBP', 'CD3D']
df_expr = pd.DataFrame(adata.raw[:, genes].X.toarray(), index=adata.obs_names, columns=genes)
df_expr['cluster'] = adata.obs['leiden_0.8'].astype(str)

mean_expr = df_expr.groupby('cluster').mean()
print("Mean expression (leiden 0.8):")
print(mean_expr[['CD3D', 'IL7R', 'CD8A', 'GNLY', 'NKG7', 'CD14', 'FCGR3A', 'MS4A1', 'FCER1A', 'PPBP']])

# Map cell types:
# Cluster 0: CD4+ T cells (CD3D+, IL7R+, CD4+)
# Cluster 1: NK / Cytotoxic T cells (GNLY+, NKG7+, CD8A+)
# Cluster 2: CD14+ Monocytes (CD14+, LYZ+, S100A9+)
# Cluster 3: B cells (MS4A1+, CD79A+)
# Cluster 4: FCGR3A+ Monocytes (FCGR3A+, LST1+, MS4A7+)
# Cluster 5: Dendritic cells (HLA-DQA1+, FCER1A+, CST3+)
# Cluster 6: Megakaryocytes / Platelets (PPBP+, PF4+)

cluster_map = {
    '0': 'CD4+ T cells',
    '1': 'NK cells',
    '2': 'CD14+ Monocytes',
    '3': 'B cells',
    '4': 'FCGR3A+ Monocytes',
    '5': 'Dendritic cells',
    '6': 'Platelets'
}

adata.obs['cell_type'] = adata.obs['leiden_0.8'].astype(str).map(cluster_map)

# Save labels.csv for filtered cells (and check if we need to include all 2700 cells or filtered cells)
labels_df = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type']
})
labels_df.to_csv('labels.csv', index=False)
print("\nSaved labels.csv with shape:", labels_df.shape)
print(labels_df['cell_type'].value_counts())
print(labels_df.head(10))