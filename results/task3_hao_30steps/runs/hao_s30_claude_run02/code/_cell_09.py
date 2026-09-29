import scanpy as sc
adata = sc.read_h5ad('adata_full_lognorm_clustered.h5ad')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
import pandas as pd
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
for g in groups:
    top = [result['names'][g][i] for i in range(10)]
    print(g, top)