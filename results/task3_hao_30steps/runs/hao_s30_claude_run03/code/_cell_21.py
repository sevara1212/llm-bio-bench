import scanpy as sc
import pandas as pd

ad = sc.read_h5ad('clustered.h5ad')

cluster_labels = {
    '0': 'CD4+ Naive T cell',
    '1': 'CD4+ Memory T cell',
    '2': 'CD8+ Effector T cell',
    '3': 'Memory B cell',
    '4': 'Proliferating NK/T cell',
    '5': 'CD8+ Effector Memory T cell',
    '6': 'cDC2 (Conventional DC)',
    '7': 'CD16+ Monocyte',
    '8': 'CD56bright NK cell',
    '9': 'Plasmacytoid Dendritic Cell (pDC)',
    '10': 'Naive B cell',
    '11': 'Plasma cell',
    '12': 'Hematopoietic Stem/Progenitor Cell',
    '13': 'CD14+ Monocyte',
    '14': 'Platelet',
    '15': 'CD56dim NK cell',
    '16': 'CD8+ Central Memory T cell',
    '17': 'cDC1 (Conventional DC)',
    '18': 'Regulatory T cell (Treg)',
    '19': 'Erythrocyte'
}

ad.obs['leiden_05'] = ad.obs['leiden_05'].astype(str)
ad.obs['cell_type'] = ad.obs['leiden_05'].map(cluster_labels)

print(ad.obs['cell_type'].value_counts())
print()
print("Any unmapped?", ad.obs['cell_type'].isna().sum())

ad.write('final_labeled.h5ad')