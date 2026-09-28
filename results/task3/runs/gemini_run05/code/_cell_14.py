# Let's inspect expression of CD4, CD8A, CD8B, GNLY, NKG7, FCGR3A, CD14 in Res 1.0 clusters 0, 1, 3
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_filtered.h5ad')
raw_df = pd.DataFrame(
    adata.raw[:, ['CD4', 'IL7R', 'CCR7', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'CD14', 'FCGR3A', 'MS4A1', 'FCER1A', 'PPBP']].X.toarray(),
    index=adata.obs_names,
    columns=['CD4', 'IL7R', 'CCR7', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'CD14', 'FCGR3A', 'MS4A1', 'FCER1A', 'PPBP']
)
raw_df['cluster_1.0'] = adata.obs['leiden_1.0']
raw_df['cluster_0.8'] = adata.obs['leiden_0.8']

print("Mean expression by cluster (res 1.0):")
print(raw_df.groupby('cluster_1.0').mean())

print("\nMean expression by cluster (res 0.8):")
print(raw_df.groupby('cluster_0.8').mean())