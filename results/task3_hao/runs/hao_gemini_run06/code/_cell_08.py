import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)

# Check resolution 0.4, 0.5, 0.6
for res in [0.3, 0.4, 0.5, 0.6]:
    sc.tl.leiden(adata, resolution=res, random_state=42, key_added=f'leiden_{res}')
    print(f"res {res}: {adata.obs[f'leiden_{res}'].nunique()} clusters")

# Let's check marker genes at resolution 0.4 or 0.5
sc.tl.rank_genes_groups(adata, 'leiden_0.4', method='wilcoxon')
marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
print("\nMarkers res 0.4:")
print(marker_genes)