# Let's inspect each cluster's top 20 marker genes and expression of key lineage markers across all clusters!
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, mask_var="highly_variable", svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42, key_added='leiden')
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

for g in sorted(groups, key=lambda x: int(x)):
    genes = [result['names'][g][i] for i in range(12)]
    print(f"Cluster {g:2s} (n={sum(adata.obs['leiden']==g):4d}): {', '.join(genes)}")