print("Data type of X:", type(adata.X))
print("Max/min value in X:", adata.X.min(), adata.X.max())

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Mitochondrial genes
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))  # Ribosomal genes

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
print(adata.obs.describe())