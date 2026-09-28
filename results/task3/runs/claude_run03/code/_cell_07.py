import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('adata_clustered.h5ad')

# Use raw (log-normalized, non-scaled) data for marker gene finding
adata_raw = adata.raw.to_adata()
adata_raw.obs['leiden'] = adata.obs['leiden']

sc.tl.rank_genes_groups(adata_raw, 'leiden', method='wilcoxon')

result = adata_raw.uns['rank_genes_groups']
groups = result['names'].dtype.names
for group in groups:
    top_genes = result['names'][group][:15]
    print(f"Cluster {group}: {list(top_genes)}")