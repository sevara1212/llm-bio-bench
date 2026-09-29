import scanpy as sc
import pandas as pd

full = sc.read_h5ad('full_with_markers.h5ad')

cluster_to_celltype = {
    '0': 'CD4+ T cell (naive/memory)',
    '1': 'CD4+ T cell (memory)',
    '2': 'CD8+ T cell (cytotoxic)',
    '3': 'B cell (memory)',
    '4': 'Proliferating lymphocyte',
    '5': 'CD4+ T cell (effector memory)',
    '6': 'Dendritic cell (cDC2)',
    '7': 'CD16+ Monocyte',
    '8': 'NK cell (CD56bright)',
    '9': 'Plasmacytoid dendritic cell (pDC)',
    '10': 'Naive B cell',
    '11': 'Plasma cell',
    '12': 'Hematopoietic stem/progenitor cell (HSPC)',
    '13': 'CD14+ Monocyte (classical)',
    '14': 'Platelet',
    '15': 'NK cell (CD56dim/cytotoxic)',
    '16': 'CD8+ T cell (memory)',
    '17': 'Dendritic cell (cDC1)',
    '18': 'Regulatory T cell (Treg)',
    '19': 'Erythrocyte',
}

full.obs['cell_type'] = full.obs['leiden'].astype(str).map(cluster_to_celltype)
print(full.obs['cell_type'].value_counts())
print(full.obs['cell_type'].isna().sum())

labels = pd.DataFrame({
    'barcode': full.obs_names,
    'cell_type': full.obs['cell_type'].values
})
labels.to_csv('labels.csv', index=False)
print(labels.head(10))
print(labels.shape)