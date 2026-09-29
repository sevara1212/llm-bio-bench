import pandas as pd
df = pd.read_csv('mean_expr_by_cluster.csv', index_col=0)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 400)
pd.set_option('display.max_rows', None)
dfT = df.T
print(dfT.loc[['GNLY','KLRD1','KLRF1','NCAM1','MS4A1']])
print()
print(dfT.loc[['LYZ','FCN1','S100A8','S100A9','VCAN','FCGR3A','MS4A7','LST1','CST3','HLA-DRA']])
print()
print(dfT.loc[['CLEC9A','CD1C','FCER1A','IL3RA','LILRA4','CLEC4C','GZMB']])