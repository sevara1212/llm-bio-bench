import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Standard filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[(adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5), :].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=40, random_state=0)

for res in [0.5, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata_hvg, resolution=res, random_state=0, key_added=f'leiden_{res}')
    print(f"Res {res}: {adata_hvg.obs[f'leiden_{res}'].nunique()} clusters")

# Let's inspect markers for res=0.6 and res=0.8
for res in [0.6, 0.8]:
    sc.tl.rank_genes_groups(adata_hvg, f'leiden_{res}', method='t-test')
    result = adata_hvg.uns['rank_genes_groups']
    print(f"\n--- Resolution {res} ---")
    for g in result['names'].dtype.names:
        top_genes = list(result['names'][g][:6])
        print(f"Cluster {g} (n={sum(adata_hvg.obs[f'leiden_{res}']==g)}): {top_genes}")