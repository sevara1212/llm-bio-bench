import scanpy as sc
adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

for q in [0.01,0.02,0.05,0.95,0.98,0.99]:
    print(q, adata.obs['n_genes_by_counts'].quantile(q), adata.obs['total_counts'].quantile(q), adata.obs['pct_counts_mt'].quantile(q))