import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('full_with_clusters.h5ad')

# Final cluster annotations based on marker gene analysis
cluster_annotations = {
    '0': 'CD4 Naive T cell',
    '1': 'CD4 Memory T cell',
    '2': 'CD8 Effector/Cytotoxic T cell',
    '3': 'Memory B cell',
    '4': 'MAIT / CD4 Effector Memory T cell',
    '5': 'CD8 Naive T cell',
    '6': 'Dendritic cell (cDC2)',
    '7': 'CD16+ Monocyte',
    '8': 'NK cell (CD56 bright)',
    '9': 'Plasmacytoid Dendritic cell (pDC)',
    '10': 'Naive B cell',
    '11': 'Plasma cell',
    '12': 'Proliferating NK/T cell',
    '13': 'Hematopoietic Progenitor cell',
    '14': 'Platelet',
    '15': 'Proliferating cell (cycling)',
    '16': 'NK cell (CD56 dim)',
    '17': 'CD14+ Monocyte',
    '18': 'CD8 Memory T cell',
    '19': 'Dendritic cell (cDC1)',
    '20': 'Activated/Regulatory-like T cell',
    '21': 'CD14+ Monocyte (activated/inflammatory)',
    '22': 'Erythrocyte',
    '23': 'Plasmacytoid Dendritic cell (pDC)',
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)

print(adata.obs['cell_type'].value_counts())
print()
print(adata.obs[['leiden','cell_type']].drop_duplicates().sort_values('leiden'))