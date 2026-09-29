# Check sample prefix distribution across clusters
batch = [b.split('_')[0] for b in adata.obs_names]
adata_proc.obs['batch'] = batch
adata_hvg.obs['batch'] = batch
print(pd.crosstab(adata_hvg.obs['batch'], adata_hvg.obs['leiden_0.5']))