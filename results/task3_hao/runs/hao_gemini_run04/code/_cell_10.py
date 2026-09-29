import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('processed.h5ad')

# Compute rank_genes_groups for leiden_0.5
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')

# Print top 5 genes for each cluster
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_genes = pd.DataFrame(
    {group: result['names'][group][:10] for group in groups}
)
print("Top genes per cluster (leiden 0.5):")
print(top_genes)