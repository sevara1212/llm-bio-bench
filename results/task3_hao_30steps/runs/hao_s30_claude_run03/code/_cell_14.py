import pandas as pd
df = pd.read_csv('mean_expr.csv', index_col=0)
for idx in df.index:
    row = df.loc[idx]
    top = row.sort_values(ascending=False).head(8)
    print(f"Cluster {idx}: {list(zip(top.index, top.round(2)))}")