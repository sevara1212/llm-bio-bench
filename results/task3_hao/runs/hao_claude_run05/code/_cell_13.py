import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('adata_final.h5ad')
adata_raw = adata.raw.to_adata()
adata_raw.obs['leiden'] = adata.obs['leiden'].values

genes_check = ['CD4','CD8A','CD8B','IL7R','CCR7','SELL','LEF1','FOXP3','IL2RA','TNFRSF18','TNFRSF4',
               'GZMK','GZMH','GZMB','GZMA','NKG7','GNLY','KLRD1','FCGR3A','KLRF1','CD3D','CD3E','TRAC',
               'MS4A1','CD79A','IGHD','TCL1A','BANK1','CD27','TCF7','CCL5','CST7']
genes_present = [g for g in genes_check if g in adata_raw.var_names]
missing = [g for g in genes_check if g not in adata_raw.var_names]
print("Missing:", missing)

df = pd.DataFrame(adata_raw[:, genes_present].X.toarray(), columns=genes_present, index=adata_raw.obs_names)
df['leiden'] = adata_raw.obs['leiden'].values
means = df.groupby('leiden')[genes_present].mean()
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)
print(means.round(2))