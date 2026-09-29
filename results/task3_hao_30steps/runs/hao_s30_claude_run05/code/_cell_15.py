import pandas as pd
df = pd.read_csv('mean_expr_by_cluster.csv', index_col=0)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 400)
pd.set_option('display.max_rows', None)
dfT = df.T
# print row by row for genes 6-30 (middle section that got truncated)
print(dfT.iloc[5:25])