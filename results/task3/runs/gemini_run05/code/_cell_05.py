# Inspect genes and duplicates
adata.var_names_make_unique()
print("Unique genes:", adata.var_names.is_unique)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(adata.obs.describe())