# Make gene names unique
adata.var_names_make_unique()

# Quality control
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']].describe())