import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('adata_clustered.h5ad')

cluster_to_celltype = {
    '0': 'CD4 T cells',
    '1': 'CD14+ Monocytes',
    '2': 'B cells',
    '3': 'CD8 T cells',
    '4': 'FCGR3A+ Monocytes',
    '5': 'NK cells',
    '6': 'Dendritic cells',
    '7': 'Megakaryocytes',
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)
print(adata.obs['cell_type'].value_counts())

labels = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type'].values
})
print(labels.head())
print(labels.shape)

labels.to_csv('labels.csv', index=False)
print("Saved labels.csv")