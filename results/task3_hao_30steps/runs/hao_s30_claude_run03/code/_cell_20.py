import scanpy as sc
import pandas as pd

ad = sc.read_h5ad('clustered.h5ad')

# Check Treg markers for cluster 18, and CD56bright vs dim markers for 8 vs 15
extra = ['FOXP3','IL2RA','TNFRSF18','TNFRSF4','CTLA4','XCL1','XCL2','KLRC1','FCGR3A','FGFBP2','GZMB','SELL']
df = sc.get.obs_df(ad, keys=extra+['leiden_05'], use_raw=True)
mean_expr = df.groupby('leiden_05').mean()
mean_expr.index = mean_expr.index.astype(str)
pd.set_option('display.width', 300)
print(mean_expr.loc[['8','15','18']].round(3))