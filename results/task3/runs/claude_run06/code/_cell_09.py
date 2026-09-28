import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('clustered.h5ad')

markers2 = ['CCR7','S100A4','SELL','IL7R','CD3D','CD8A','CD8B','GZMK','CD27']
markers2 = [m for m in markers2 if m in adata.raw.var_names]
df = pd.DataFrame(adata.raw[:, markers2].X.toarray(), columns=markers2, index=adata.obs_names)
df['leiden'] = adata.obs['leiden'].values
print(df.groupby('leiden').mean().round(2))

print(adata.obs['leiden'].value_counts())