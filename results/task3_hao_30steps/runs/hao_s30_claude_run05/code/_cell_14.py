import pandas as pd
df = pd.read_csv('mean_expr_by_cluster.csv', index_col=0)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 400)
print(df.T)  # transposed to see full - genes as rows, clusters as columns