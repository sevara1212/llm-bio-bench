# Check mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.str.startswith(('HBA', 'HBB', 'HBD', 'HBE', 'HBG', 'HBM', 'HBQ', 'HBZ'))

print(f"Mitochondrial genes count: {adata.var['mt'].sum()}")
print(f"Ribosomal genes count: {adata.var['ribo'].sum()}")
print(f"Hemoglobin genes count: {adata.var['hb'].sum()}")

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, log1p=False, inplace=True)
print(adata.obs.describe())