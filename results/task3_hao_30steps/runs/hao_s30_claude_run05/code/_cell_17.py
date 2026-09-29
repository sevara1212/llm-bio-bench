import pandas as pd
df = pd.read_csv('mean_expr_by_cluster.csv', index_col=0)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 400)
dfT = df.T
for g in ['LYZ','FCN1','S100A8','S100A9','VCAN','FCGR3A','MS4A7','LST1']:
    print(g)
    print(dfT.loc[g])
    print()