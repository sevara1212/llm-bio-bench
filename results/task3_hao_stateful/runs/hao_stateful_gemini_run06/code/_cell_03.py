# Let's inspect QC metrics
# mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
# ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
# hemoglobin genes
adata.var['hb'] = adata.var_names.str.startswith(('HBA', 'HBB', 'HBD', 'HBE', 'HBG', 'HBM', 'HBQ', 'HBZ'))

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, log1p=False, inplace=True)
print(adata.obs.describe())