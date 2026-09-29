import scanpy as sc
import pandas as pd

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
sc.tl.leiden(adata, resolution=0.5, random_state=42)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

for g in groups:
    names = [result['names'][g][i] for i in range(10)]
    print(f"Cluster {g} (n={sum(adata.obs['leiden']==g)}): {', '.join(names)}")