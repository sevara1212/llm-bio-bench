import scanpy as sc
adata = sc.read_h5ad('adata_processed.h5ad')
adata_hvg = sc.read_h5ad('adata_hvg.h5ad')

adata.obs['leiden'] = adata_hvg.obs['leiden_0.4'].values
print(adata.obs['leiden'].value_counts())

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', use_raw=True)
import pandas as pd
res = adata.uns['rank_genes_groups']
groups = res['names'].dtype.names
top_markers = {}
for g in groups:
    top_markers[g] = [res['names'][g][i] for i in range(15)]
df = pd.DataFrame(top_markers)
print(df)
df.to_csv('top_markers.csv')
adata.write('adata_processed.h5ad')