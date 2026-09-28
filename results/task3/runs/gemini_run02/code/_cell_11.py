import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_processed.h5ad')
sc.tl.leiden(adata, resolution=1.0, key_added='leiden')
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_genes = pd.DataFrame(
    {group: result['names'][group][:10] for group in groups}
)
print("Top 10 markers per cluster (res=1.0):")
print(top_genes)