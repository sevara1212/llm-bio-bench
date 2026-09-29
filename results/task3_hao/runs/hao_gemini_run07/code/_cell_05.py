import scanpy as sc
adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print(adata.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']].describe(percentiles=[0.005, 0.01, 0.05, 0.5, 0.95, 0.98, 0.99, 0.995]))