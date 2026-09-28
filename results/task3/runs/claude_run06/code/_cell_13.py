import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('clustered.h5ad')

cluster_to_celltype = {
    '0': 'CD4 T cells',
    '1': 'CD14+ Monocytes',
    '2': 'B cells',
    '3': 'CD8 T cells',
    '4': 'FCGR3A+ Monocytes',
    '5': 'NK cells',
    '6': 'CD4 T cells',
    '7': 'Dendritic cells',
    '8': 'Megakaryocytes',
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)

labels = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type'].values
})

labels.to_csv('labels.csv', index=False)
print(labels.shape)
print(labels.head(10))
print(labels['cell_type'].value_counts())