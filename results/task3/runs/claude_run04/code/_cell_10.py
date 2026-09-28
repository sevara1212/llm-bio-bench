import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('clustered.h5ad')

cluster_to_celltype = {
    '0': 'CD4+ T cells',
    '1': 'B cells',
    '2': 'CD16+ Monocytes',
    '3': 'CD14+ Monocytes',
    '4': 'NK cells',
    '5': 'CD8+ T cells',
    '6': 'Dendritic cells',
    '7': 'Megakaryocytes',
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)
print(adata.obs['cell_type'].value_counts())

labels = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type'].values
})
labels.to_csv('labels.csv', index=False)
print(labels.head())
print(labels.shape)