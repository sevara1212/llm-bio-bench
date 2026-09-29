import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('adata_full_annotated.h5ad')

cluster_to_celltype = {
    '0': 'CD4+ Naive T cell',
    '1': 'CD4+ Memory T cell',
    '2': 'CD8+ Effector/Cytotoxic T cell',
    '3': 'Memory B cell',
    '4': 'Proliferating T/NK cell',
    '5': 'CD8+ Effector Memory T cell',
    '6': 'Dendritic cell (cDC2)',
    '7': 'CD16+ Monocyte',
    '8': 'NK cell (CD56bright)',
    '9': 'Plasmacytoid dendritic cell (pDC)',
    '10': 'Naive B cell',
    '11': 'Plasma cell',
    '12': 'Hematopoietic progenitor cell',
    '13': 'CD14+ Monocyte',
    '14': 'Platelet',
    '15': 'NK cell (CD56dim/cytotoxic)',
    '16': 'CD8+ Memory T cell',
    '17': 'Dendritic cell (cDC1)',
    '18': 'Regulatory T cell (Treg)',
    '19': 'Erythrocyte (contaminant)'
}

adata.obs['cell_type'] = adata.obs['leiden_r05'].astype(str).map(cluster_to_celltype)
print(adata.obs['cell_type'].value_counts())
print(adata.obs['cell_type'].isna().sum())