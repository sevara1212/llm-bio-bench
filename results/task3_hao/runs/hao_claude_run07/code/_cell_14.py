import pandas as pd
mean_expr = pd.read_csv('marker_mean_expr.csv', index_col=0)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 250)
print(mean_expr.round(2))