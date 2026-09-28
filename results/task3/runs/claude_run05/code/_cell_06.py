import scanpy as sc
adata = sc.read_h5ad('adata_filtered.h5ad')

# save raw counts
adata.layers['counts'] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# save normalized/log data for marker gene analysis later
adata.raw = adata

# HVG
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Number of HVGs:", adata.var.highly_variable.sum())

adata = adata[:, adata.var.highly_variable].copy()

# regress out and scale
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

adata.write('adata_hvg_scaled.h5ad')
print(adata)