import scanpy as sc
adata = sc.read_h5ad('clustered.h5ad')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

import pandas as pd
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

for grp in groups:
    top_genes = [result['names'][grp][i] for i in range(15)]
    print(f"Cluster {grp}: {top_genes}")