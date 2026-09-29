import pandas as pd
pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 50)
df = pd.read_csv('mean_expr.csv', index_col=0)
print(df.T)