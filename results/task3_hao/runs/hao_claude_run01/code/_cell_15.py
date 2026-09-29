import scanpy as sc
import pandas as pd

ad = sc.read_h5ad('processed.h5ad')

cluster_to_celltype = {
    '0': 'CD8+ Effector Memory T cell',
    '1': 'CD4+ Memory T cell',
    '2': 'CD4+ Naive T cell',
    '3': 'Memory B cell',
    '4': 'Proliferating Lymphocyte',
    '5': 'CD8+ Naive T cell',
    '6': 'MAIT/GZMK+ T cell',
    '7': 'CD16+ Monocyte',
    '8': 'Conventional Dendritic Cell (cDC2)',
    '9': 'Naive B cell',
    '10': 'NK cell (CD56bright)',
    '11': 'Plasmacytoid Dendritic Cell (pDC)',
    '12': 'Hematopoietic Stem/Progenitor Cell',
    '13': 'Platelet',
    '14': 'Plasma cell',
    '15': 'NK cell (CD56dim)',
    '16': 'CD8+ Memory T cell (GZMK+)',
    '17': 'CD14+ Monocyte (Classical)',
    '18': 'Conventional Dendritic Cell (cDC1)',
    '19': 'Regulatory T cell (Treg)',
    '20': 'CD14+ Monocyte (Inflammatory)',
    '21': 'Plasmacytoid Dendritic Cell (pDC)',
    '22': 'Erythrocyte',
    '23': 'Erythrocyte',
}

ad.obs['cell_type'] = ad.obs['leiden'].astype(str).map(cluster_to_celltype)
print(ad.obs['cell_type'].value_counts())
print(ad.obs['cell_type'].isna().sum())

labels = ad.obs[['cell_type']].copy()
labels.index.name = 'barcode'
labels = labels.reset_index()
labels.to_csv('labels.csv', index=False)
print(labels.head(10))
print(labels.shape)