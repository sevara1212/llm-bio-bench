import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('adata_final.h5ad')

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)

top_markers = {}
for group in groups:
    top_markers[group] = [result['names'][group][i] for i in range(20)]

df_markers = pd.DataFrame(top_markers)
print(df_markers.to_string())