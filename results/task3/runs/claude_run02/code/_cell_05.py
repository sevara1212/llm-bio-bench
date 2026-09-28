import scanpy as sc
adata = sc.read_h5ad('adata_filtered.h5ad')

# Save raw counts
adata.layers['counts'] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(adata.var.highly_variable.sum())

adata.raw = adata  # store full data for later marker gene lookup

adata_hvg = adata[:, adata.var.highly_variable].copy()

# Regress out and scale
sc.pp.scale(adata_hvg, max_value=10)

adata_hvg.write('adata_hvg.h5ad')
adata.write('adata_norm.h5ad')