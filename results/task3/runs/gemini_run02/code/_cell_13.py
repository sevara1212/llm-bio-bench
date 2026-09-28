import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_processed.h5ad')
sc.tl.leiden(adata, resolution=1.0, key_added='leiden')

markers = ['IL7R', 'CD3D', 'CD4', 'CD8A', 'CD8B', 'S100A8', 'CD14', 'FCGR3A', 'MS4A7', 
           'MS4A1', 'CD79A', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP']

mean_df = pd.DataFrame(index=markers, columns=[str(i) for i in range(9)])
for g in markers:
    if g in adata.raw.var_names:
        expr = adata.raw[:, g].X.toarray().flatten()
        for c in range(9):
            mean_df.loc[g, str(c)] = np.round(expr[adata.obs['leiden'] == str(c)].mean(), 3)

print("Marker mean expression per cluster:")
print(mean_df)
print("\nCluster counts:")
print(adata.obs['leiden'].value_counts().sort_index())