import scanpy as sc
adata = sc.read_h5ad('qc.h5ad')

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", adata.var.highly_variable.sum())

adata_hvg = adata[:, adata.var.highly_variable].copy()

sc.pp.regress_out(adata_hvg, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata_hvg, max_value=10)

adata_hvg.write('processed.h5ad')
print(adata_hvg)