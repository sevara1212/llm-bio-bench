# Let's inspect raw.to_adata() or adata.raw[:, gene].X
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

# Let's check Leiden at resolution 1.0 vs 0.8 vs 0.5
adata_sub = adata.copy()
sc.pp.highly_variable_genes(adata_sub, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_sub = adata_sub[:, adata_sub.var.highly_variable].copy()
sc.pp.scale(adata_sub, max_value=10)
sc.pp.pca(adata_sub, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_sub, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.leiden(adata_sub, resolution=1.0, key_added='leiden_1', random_state=0)

adata.obs['leiden_1'] = adata_sub.obs['leiden_1']
sc.tl.rank_genes_groups(adata, 'leiden_1', method='wilcoxon')
result = adata.uns['rank_genes_groups']
print("Markers res 1.0:")
for g in sorted(result['names'].dtype.names, key=lambda x: int(x)):
    top_genes = list(result['names'][g][:8])
    print(f"Cluster {g} (n={sum(adata.obs['leiden_1']==g)}): {top_genes}")