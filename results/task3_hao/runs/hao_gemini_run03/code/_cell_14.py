import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('processed.h5ad')
adata.raw = None
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon', use_raw=False)

for i in range(21):
    top_diff = list(adata.uns['rank_genes_groups']['names'][str(i)][:7])
    print(f"Cluster {i} (n={sum(adata.obs['leiden_0.5']==str(i))}): {top_diff}")