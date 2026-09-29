import pandas as pd
mean_expr = pd.read_csv('marker_mean_expr.csv', index_col=0)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 250)
print(mean_expr.loc[['MS4A1','CD79A','CD79B','BANK1','CD14','LYZ','S100A8','S100A9','FCN1','FCGR3A']].round(2))
print()
print(adata_cluster_sizes if False else "")