import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_filtered.h5ad')

# Check markers for res 0.6 / 0.8
for res in ['leiden_0.6', 'leiden_0.8']:
    print(f"=== Results for {res} ===")
    sc.tl.rank_genes_groups(adata, groupby=res, method='wilcoxon')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df = pd.DataFrame({group: [f"{result['names'][group][i]} ({result['scores'][group][i]:.1f})" for i in range(8)] for group in groups})
    print(df)