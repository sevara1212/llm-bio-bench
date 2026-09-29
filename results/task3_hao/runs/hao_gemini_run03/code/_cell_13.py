import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed.h5ad')

for i in range(21):
    c_adata = adata[adata.obs['leiden_0.5'] == str(i)]
    # calculate top expressing genes or differential genes vs rest
    top_diff = adata.uns['rank_genes_groups']['names'][str(i)][:5]
    print(f"Cluster {i} (n={c_adata.n_obs}): top markers = {list(top_diff)}")