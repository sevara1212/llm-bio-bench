import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_sub = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_sub, max_value=10)
sc.tl.pca(adata_sub, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_sub, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.umap(adata_sub, random_state=0)

# Let's test a few resolutions: 0.5, 0.8, 1.0
for res in [0.4, 0.5, 0.6, 0.8]:
    sc.tl.leiden(adata_sub, resolution=res, key_added=f'leiden_{res}', random_state=0)
    print(f"Res {res}: {adata_sub.obs[f'leiden_{res}'].nunique()} clusters")

adata.obs['leiden'] = adata_sub.obs['leiden_0.5']
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
result = adata.uns['rank_genes_groups']
for group in sorted(result['names'].dtype.names, key=lambda x: int(x)):
    print(f"\n--- Cluster {group} (n={sum(adata.obs['leiden'] == group)}) ---")
    top_genes = result['names'][group][:15]
    top_scores = result['scores'][group][:15]
    print(list(zip(top_genes, [round(s, 2) for s in top_scores])))