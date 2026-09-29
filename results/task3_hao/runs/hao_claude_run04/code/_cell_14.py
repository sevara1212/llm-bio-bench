import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('clustered.h5ad')

cluster_to_celltype = {
    '0': 'Naive CD4 T cells',
    '1': 'Memory CD4 T cells',
    '2': 'CD8+ Cytotoxic T cells',
    '3': 'Memory B cells',
    '4': 'Proliferating lymphocytes',
    '5': 'Conventional dendritic cells (cDC2)',
    '6': 'MAIT cells',
    '7': 'CD16+ NK cells',
    '8': 'CD16+ Monocytes',
    '9': 'CD56bright NK cells',
    '10': 'Naive B cells',
    '11': 'Plasmacytoid dendritic cells',
    '12': 'Hematopoietic stem/progenitor cells',
    '13': 'CD14+ Monocytes',
    '14': 'Platelets',
    '15': 'Plasma cells',
    '16': 'Conventional dendritic cells (cDC1)',
    '17': 'Memory CD8 T cells',
    '18': 'Gamma-delta T cells',
    '19': 'Erythrocytes',
}

adata.obs['cell_type'] = adata.obs['leiden_05'].astype(str).map(cluster_to_celltype)
print(adata.obs['cell_type'].value_counts())
print(adata.obs['cell_type'].isna().sum())

labels = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type'].values
})
labels.to_csv('labels.csv', index=False)
print(labels.head())
print(labels.shape)