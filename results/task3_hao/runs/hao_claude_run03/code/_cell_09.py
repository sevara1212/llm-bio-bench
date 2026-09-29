import scanpy as sc
import pandas as pd

adata_hvg = sc.read_h5ad('clustered.h5ad')

# Build full adata with raw log-normalized data + leiden labels, for marker gene analysis on ALL genes
full = adata_hvg.raw.to_adata()
full.obs['leiden'] = adata_hvg.obs['leiden'].values

sc.tl.rank_genes_groups(full, 'leiden', method='wilcoxon')

full.write('full_with_markers.h5ad')

# print top markers per cluster
result = full.uns['rank_genes_groups']
groups = result['names'].dtype.names
for g in groups:
    top = [result['names'][g][i] for i in range(15)]
    print(g, top)