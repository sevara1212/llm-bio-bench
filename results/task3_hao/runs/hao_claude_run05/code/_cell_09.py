import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('adata_clustered.h5ad')

# Use leiden_05 as main clustering
adata.obs['leiden'] = adata.obs['leiden_05']

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', use_raw=True)

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

top_markers = {}
for group in groups:
    top_markers[group] = [result['names'][group][i] for i in range(15)]

df_markers = pd.DataFrame(top_markers)
print(df_markers)

adata.write('adata_final.h5ad')