import scanpy as sc
import pandas as pd

ad = sc.read_h5ad('clustered.h5ad')
extra = ['CD3D','CD3E','CD4','CD8A','CD8B','IL7R','CCR7','SELL','GZMK','GZMB','NKG7','GNLY',
         'KLRD1','KLRF1','FCGR3A','KLRB1','CD27','GPR183','TCF7','MKI67','STMN1','FOXP3','IL2RA']

df = sc.get.obs_df(ad, keys=extra+['leiden_05'], use_raw=True)
mean_expr = df.groupby('leiden_05').mean()
mean_expr.index = mean_expr.index.astype(str)

tcell_clusters = ['0','1','2','4','5','8','15','16','18']
pd.set_option('display.width', 300)
print(mean_expr.loc[tcell_clusters].round(2))