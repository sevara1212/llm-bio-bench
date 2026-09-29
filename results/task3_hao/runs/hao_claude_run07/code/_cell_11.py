import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('adata_processed.h5ad')

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', use_raw=True)
res = adata.uns['rank_genes_groups']
groups = res['names'].dtype.names
top_markers = {}
for g in groups:
    top_markers[g] = [res['names'][g][i] for i in range(15)]
df = pd.DataFrame(top_markers)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
print(df)
df.to_csv('top_markers.csv')
adata.write('adata_processed.h5ad')