import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('processed_adata.h5ad')

# Check expression of key PBMC markers across leiden_0.9 and leiden_1.0
genes = ['CD3D', 'IL7R', 'CD4', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'MS4A1', 'CD14', 'FCGR3A', 'CST3', 'PPBP', 'FCER1A']

# Mean expression per cluster in leiden_0.9:
raw_df = pd.DataFrame(adata.raw[:, genes].X.toarray(), columns=genes, index=adata.obs_names)
raw_df['leiden_0.9'] = adata.obs['leiden_0.9'].values

mean_exp = raw_df.groupby('leiden_0.9').mean()
print("Mean expression per cluster in leiden_0.9:")
print(mean_exp.round(2))