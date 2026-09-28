import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_sub = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_sub, max_value=10)
sc.pp.pca(adata_sub, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_sub, n_neighbors=10, n_pcs=40, random_state=0)

# Try different leiden resolutions
for res in [0.4, 0.5, 0.6, 0.8]:
    sc.tl.leiden(adata_sub, resolution=res, key_added=f'leiden_{res}', random_state=0)
    print(f"Resolution {res}: {adata_sub.obs[f'leiden_{res}'].nunique()} clusters")

# Let's check markers for resolution 0.5 or 0.6
sc.tl.leiden(adata_sub, resolution=0.5, key_added='leiden', random_state=0)
adata.obs['leiden'] = adata_sub.obs['leiden']
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markers = pd.DataFrame({group: result['names'][group][:10] for group in groups})
print("\nTop 10 markers per cluster (res=0.5):")
print(markers)