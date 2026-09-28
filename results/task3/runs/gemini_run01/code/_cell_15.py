import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed_adata.h5ad')

# Mapping for leiden_0.9:
cluster_to_celltype = {
    '0': 'Memory CD4+ T cell',
    '1': 'Naive CD4+ T cell',
    '2': 'CD14+ Monocyte',
    '3': 'B cell',
    '4': 'CD8+ T cell',
    '5': 'FCGR3A+ Monocyte',
    '6': 'NK cell',
    '7': 'Dendritic cell',
    '8': 'Platelet'
}

adata.obs['cell_type'] = adata.obs['leiden_0.9'].map(cluster_to_celltype)

df_labels = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type'].values
})

df_labels.to_csv('labels.csv', index=False)

print("Saved labels.csv with shape:", df_labels.shape)
print("\nFirst 10 rows:")
print(df_labels.head(10))
print("\nCell type distribution:")
print(df_labels['cell_type'].value_counts())