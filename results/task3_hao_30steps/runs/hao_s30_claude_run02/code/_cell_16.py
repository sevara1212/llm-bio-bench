import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('adata_full_lognorm_clustered.h5ad')

cluster_to_celltype = {
    '0': 'Naive CD4+ T cells',
    '1': 'Memory CD4+ T cells',
    '2': 'CD8+ effector/cytotoxic T cells',
    '3': 'Memory B cells',
    '4': 'Proliferating lymphocytes',
    '5': 'MAIT cells',
    '6': 'Conventional dendritic cells (cDC2)',
    '7': 'CD16+ non-classical monocytes',
    '8': 'CD56bright NK cells',
    '9': 'Plasmacytoid dendritic cells (pDC)',
    '10': 'Naive B cells',
    '11': 'Plasma cells',
    '12': 'Hematopoietic stem/progenitor cells (HSPC)',
    '13': 'CD14+ classical monocytes',
    '14': 'Platelets/Megakaryocytes',
    '15': 'CD56dim NK cells',
    '16': 'Memory CD8+ T cells',
    '17': 'Conventional dendritic cells (cDC1)',
    '18': 'Regulatory T cells (Tregs)',
    '19': 'Erythroid cells',
}

adata.obs['cell_type'] = adata.obs['leiden'].astype(str).map(cluster_to_celltype)
print(adata.obs['cell_type'].value_counts())
print(adata.obs['cell_type'].isna().sum())

labels = pd.DataFrame({'barcode': adata.obs_names, 'cell_type': adata.obs['cell_type'].values})
labels.to_csv('labels.csv', index=False)
print(labels.head())
print(labels.shape)