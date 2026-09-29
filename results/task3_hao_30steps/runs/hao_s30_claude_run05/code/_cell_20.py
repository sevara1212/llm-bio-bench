import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('full_with_clusters.h5ad')

cluster_to_celltype = {
    '0': 'Naive CD4+ T cell',
    '1': 'Memory CD4+ T cell',
    '2': 'CD8+ Cytotoxic T cell',
    '3': 'Memory B cell',
    '4': 'CD4+ Effector Memory T cell',
    '5': 'Naive CD8+ T cell',
    '6': 'Conventional Dendritic Cell (cDC2)',
    '7': 'CD16+ Monocyte (Non-classical)',
    '8': 'NK cell (CD56bright)',
    '9': 'Plasmacytoid Dendritic Cell (pDC)',
    '10': 'Naive B cell',
    '11': 'Plasma cell',
    '12': 'Proliferating NK cell',
    '13': 'Hematopoietic Stem/Progenitor Cell (HSPC)',
    '14': 'Platelet/Megakaryocyte',
    '15': 'Proliferating T cell',
    '16': 'NK cell (CD56dim, cytotoxic)',
    '17': 'Classical Monocyte (CD14+)',
    '18': 'CD8+ Memory T cell',
    '19': 'Conventional Dendritic Cell (cDC1)',
    '20': 'Regulatory T cell (Treg)',
    '21': 'Classical Monocyte (activated/inflammatory)',
    '22': 'Erythroid cell',
    '23': 'Plasmacytoid Dendritic Cell (pDC)',
}

adata.obs['cell_type'] = adata.obs['leiden'].astype(str).map(cluster_to_celltype)
print(adata.obs['cell_type'].value_counts())
print(adata.obs['cell_type'].isna().sum())