import pandas as pd
pd.set_option('display.width', 300)
pd.set_option('display.max_columns', 60)
df = pd.read_csv('mean_expr.csv', index_col=0)
print(df.to_string())