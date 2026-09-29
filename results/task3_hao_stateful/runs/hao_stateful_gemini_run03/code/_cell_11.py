# Wait, why are there so many clusters? Could batch effect across samples (L1-L5, E2L1-E2L8) be causing clusters?
# Let's check distribution of prefixes across clusters!
adata_proc.obs['sample'] = adata_proc.obs_names.map(lambda x: x.split('_')[0])
print(pd.crosstab(adata_proc.obs['leiden_0.5'], adata_proc.obs['sample']))