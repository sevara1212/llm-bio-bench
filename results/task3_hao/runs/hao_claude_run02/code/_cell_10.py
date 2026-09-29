import scanpy as sc
adata = sc.read_h5ad('clustered2.h5ad')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=False)

import pandas as pd
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_markers = {}
for g in groups:
    top_markers[g] = [result['names'][g][i] for i in range(15)]

for g in groups:
    print(g, top_markers[g])