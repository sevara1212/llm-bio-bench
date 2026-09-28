import scanpy as sc
import pandas as pd

adata_hvg = sc.read_h5ad('clustered.h5ad')
adata_full = sc.read_h5ad('processed_full.h5ad')

# transfer cluster labels to full (log-normalized, unscaled) data for marker gene testing
adata_full.obs['leiden'] = adata_hvg.obs['leiden'].values

sc.tl.rank_genes_groups(adata_full, 'leiden', method='wilcoxon')

for cl in adata_full.obs['leiden'].cat.categories:
    names = adata_full.uns['rank_genes_groups']['names'][cl][:15]
    print(cl, list(names))