import pandas as pd
df = pd.read_csv('mean_expr_by_cluster.csv', index_col=0)
dfT = df.T
genes = ['LYZ','FCN1','S100A8','S100A9','VCAN','FCGR3A','MS4A7']
sub = dfT.loc[genes]
sub.to_csv('sub_check.csv')
print(sub.to_string())