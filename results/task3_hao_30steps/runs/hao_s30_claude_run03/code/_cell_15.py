import pandas as pd
df = pd.read_csv('mean_expr.csv', index_col=0)
cols = ['CD3D','CD3E','CD4','CD8A','CD8B','IL7R','CCR7','SELL','GZMK','GZMB','GZMH','NKG7','GNLY','KLRD1','KLRF1','FCGR3A']
print(df.loc[[0,1,5,16,18], cols])