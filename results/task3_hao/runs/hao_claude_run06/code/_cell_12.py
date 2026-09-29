import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('clustered_final.h5ad')

sc.tl.rank_genes_groups(adata, 'leiden_final', method='wilcoxon')

for cl in adata.obs['leiden_final'].cat.categories:
    names = adata.uns['rank_genes_groups']['names'][cl][:15]
    print(cl, list(names))