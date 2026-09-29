import pandas as pd
m=pd.read_csv('cluster_markers_top30.csv')
for c in m:
 print(c,':',', '.join(m[c].head(20).tolist()))