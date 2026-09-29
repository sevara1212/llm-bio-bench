import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('adata_full_lognorm_clustered.h5ad')

markers = ['CD34','PRSS57','ANKRD28','FOXP3','IL2RA','CTLA4','TNFRSF18','TNFRSF4','TRDC','TRGC1','TRGC2','KLRC1','XCL1','XCL2','CD27','GPR183','TRAC']
markers = [m for m in markers if m in adata.var_names]
df = sc.get.obs_df(adata, keys=markers+['leiden'])
mean_expr = df.groupby('leiden').mean()
pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 60)
print(mean_expr.round(2).T)