import scanpy as sc
import pandas as pd

adata_hvg = sc.read_h5ad('clustered.h5ad')

cluster_to_celltype = {
    '0': 'CD4+ Memory T cells',
    '1': 'CD8+ T cells',
    '2': 'NK cells',
    '3': 'CD4+ Naive T cells',
    '4': 'B cells',
    '5': 'FCGR3A+ Monocytes',
    '6': 'CD14+ Monocytes',
    '7': 'Dendritic cells',
    '8': 'Megakaryocytes/Platelets',
}

adata_hvg.obs['cell_type'] = adata_hvg.obs['leiden'].map(cluster_to_celltype)

print(adata_hvg.obs['cell_type'].value_counts())
print(adata_hvg.obs[['leiden','cell_type']].head())

labels = pd.DataFrame({
    'barcode': adata_hvg.obs_names,
    'cell_type': adata_hvg.obs['cell_type'].values
})
labels.to_csv('labels.csv', index=False)
print(labels.head())
print(labels.shape)