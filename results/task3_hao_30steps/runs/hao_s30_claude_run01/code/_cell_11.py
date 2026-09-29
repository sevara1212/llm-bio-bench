import scanpy as sc

adata = sc.read_h5ad('full_with_clusters.h5ad')

# Rank genes for clusters using log-normalized full data (not scaled)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# print top markers for each cluster
import pandas as pd
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
for group in groups:
    top_genes = result['names'][group][:15]
    print(f"Cluster {group}: {list(top_genes)}")