import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed.h5ad')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

for g in groups:
    genes = result['names'][g][:8]
    scores = result['scores'][g][:8]
    print(f"Cluster {g} (n={sum(adata.obs['leiden_0.5'] == g)}): {', '.join([f'{gene}' for gene in genes])}")