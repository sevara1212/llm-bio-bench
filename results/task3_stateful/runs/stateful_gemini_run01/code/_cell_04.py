adata.var_names_make_unique()
print(adata)
# Check mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
print(f"Number of MT genes: {adata.var['mt'].sum()}")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print(adata.obs.describe())