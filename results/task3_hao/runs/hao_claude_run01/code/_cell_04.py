import scanpy as sc
ad = sc.read_h5ad('raw_counts.h5ad')
ad.var['mt'] = ad.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

for q in [0.01,0.02,0.05,0.95,0.98,0.99]:
    print(q, ad.obs['n_genes_by_counts'].quantile(q), ad.obs['total_counts'].quantile(q), ad.obs['pct_counts_mt'].quantile(q))