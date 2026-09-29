# Let's inspect cluster 0, 1, 9, 11 more closely to be 100% sure of cell types
# Also write labels.csv
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=3)

# Log-normalized layer
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata, resolution=0.3, key_added='leiden', random_state=42)

# Check clusters 0 and 1:
for cl in ['0', '1', '9', '11']:
    print(f"\nCluster {cl}:")
    for gene in ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'IL7R', 'CCR7', 'S100A9', 'CD14', 'FCER1A', 'CLEC9A', 'WDFY4', 'MPO', 'ELANE', 'GATA2']:
        if gene in adata.raw.var_names:
            val = np.expm1(adata.raw[adata.obs['leiden'] == cl, gene].X.toarray()).mean()
            print(f"  {gene}: {val:.3f}")