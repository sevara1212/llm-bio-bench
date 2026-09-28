# Let's inspect standard scanpy tutorial on PBMC3k
# Standard tutorial:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

# Let's run this pipeline and see clusters and marker genes
import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
raw_cells = adata.obs_names.copy()

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

print(f"Remaining cells: {adata.n_obs}")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_hvg = adata[:, adata.var.highly_variable]
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.leiden(adata_hvg, resolution=0.5, random_state=0)

adata.obs['leiden'] = adata_hvg.obs['leiden']

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markers_df = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']}
)
print("Clusters:", adata.obs['leiden'].value_counts())
for g in groups:
    top_genes = list(result['names'][g][:8])
    print(f"Cluster {g}: {top_genes}")