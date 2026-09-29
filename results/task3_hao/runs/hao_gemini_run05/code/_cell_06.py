import scanpy as sc

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print(adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']].quantile([0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99]))
print("Cells with pct_counts_mt > 15:", (adata.obs['pct_counts_mt'] > 15).sum())
print("Cells with pct_counts_mt > 20:", (adata.obs['pct_counts_mt'] > 20).sum())
print("Cells with n_genes < 200:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("Cells with n_genes < 500:", (adata.obs['n_genes_by_counts'] < 500).sum())
print("Min n_genes:", adata.obs['n_genes_by_counts'].min())
print("Max n_genes:", adata.obs['n_genes_by_counts'].max())