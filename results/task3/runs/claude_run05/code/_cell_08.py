import scanpy as sc
adata = sc.read_h5ad('adata_clustered.h5ad')

# Use raw (normalized, log1p, non-scaled) data for marker gene detection
adata_full = adata.raw.to_adata()
adata_full.obs['leiden'] = adata.obs['leiden']

sc.tl.rank_genes_groups(adata_full, 'leiden', method='wilcoxon')

import pandas as pd
result = adata_full.uns['rank_genes_groups']
groups = result['names'].dtype.names
for group in groups:
    top_genes = result['names'][group][:15]
    print(f"Cluster {group}: {list(top_genes)}")