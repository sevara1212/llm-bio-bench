# Let's inspect the 8 clusters in scanpy pbmc3k tutorial:
# In scanpy tutorial:
# marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
#                 'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
#                 'FCGR3A', 'MS4A7', 'CST3', 'FCER1A', 'PPBP']
# Cell types in PBMC 3k:
# ['CD4 T', 'CD14+ Monocytes', 'B', 'CD8 T', 'NK', 'FCGR3A+ Monocytes', 'Dendritic', 'Megakaryocytes']

adata.obs['leiden'] = adata.obs['leiden_1.0']
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
import pandas as pd
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_genes = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'logfoldchanges', 'pvals_adj']}).head(10)
print(top_genes[[c for c in top_genes.columns if '_n' in c]])