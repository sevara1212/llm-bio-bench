# Let's inspect mitochondrial genes, ribosomal genes, QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
print("Number of MT genes:", adata.var['mt'].sum())
print("Number of Ribo genes:", adata.var['ribo'].sum())

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
print(adata.obs.describe())