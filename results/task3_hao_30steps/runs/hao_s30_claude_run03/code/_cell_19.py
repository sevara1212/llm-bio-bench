import pandas as pd
pd.set_option('display.width', 300)
pd.set_option('display.max_columns', 30)
df = pd.read_csv('mean_expr.csv', index_col=0)  # doesn't have all cols though

# Let's just print full mean_expr.csv with all columns for T-cell subclusters
print(df.columns.tolist())