import scanpy as sc
adata_hvg = sc.read_h5ad('adata_clustered.h5ad')
adata = sc.read_h5ad('adata_norm.h5ad')

# transfer leiden labels
adata.obs['leiden'] = adata_hvg.obs['leiden'].values

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

import pandas as pd
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
for g in groups:
    top_genes = [result['names'][g][i] for i in range(15)]
    print(f"Cluster {g}: {top_genes}")

adata.write('adata_ranked.h5ad')