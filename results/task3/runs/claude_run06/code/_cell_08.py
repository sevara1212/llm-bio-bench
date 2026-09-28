import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('clustered.h5ad')

markers = ['CD3D','CD3E','CD4','CD8A','IL7R','CCR7','CD14','LYZ','FCGR3A','MS4A7','FCER1A','CST3',
           'MS4A1','CD79A','GNLY','NKG7','PPBP','PF4']

markers = [m for m in markers if m in adata.raw.var_names]
df = pd.DataFrame(adata.raw[:, markers].X.toarray(), columns=markers, index=adata.obs_names)
df['leiden'] = adata.obs['leiden'].values

mean_expr = df.groupby('leiden').mean()
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 30)
print(mean_expr.round(2))