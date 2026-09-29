# Ah, adata[adata.obs['leiden'] == c, g].raw.X was taking all genes! In AnnData, .raw.X is (n_obs, n_raw_vars)
# Let's write a proper check and assign labels
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('processed.h5ad')

# Let's inspect the top markers for each cluster thoroughly
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
markers_df = pd.DataFrame()
for i in adata.obs['leiden'].cat.categories:
    markers_df[f"Cluster_{i}"] = adata.uns['rank_genes_groups']['names'][i][:15]
print("Top 15 markers per cluster:")
print(markers_df.T.to_string())