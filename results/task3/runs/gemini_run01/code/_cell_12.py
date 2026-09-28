import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed_adata.h5ad')

# In Scanpy's official PBMC tutorial:
# marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
#                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
#                'FCGR3A', 'MS4A7', 'CST3', 'PPBP']
# tutorial labels:
# ['CD4 T', 'CD14+ Monocytes', 'B', 'CD8 T', 'NK', 'FCGR3A+ Monocytes', 'Dendritic', 'Megakaryocytes']
# (Note: at res 0.5 there were 8 clusters in older scanpy, or let's check what clusters were found)

# Let's check diff between cluster 0 and cluster 1 directly:
sc.tl.rank_genes_groups(adata, 'leiden_0.9', groups=['0'], reference='1', method='wilcoxon', use_raw=True)
print("0 vs 1:", adata.uns['rank_genes_groups']['names']['0'][:10])

sc.tl.rank_genes_groups(adata, 'leiden_0.9', groups=['1'], reference='0', method='wilcoxon', use_raw=True)
print("1 vs 0:", adata.uns['rank_genes_groups']['names']['1'][:10])