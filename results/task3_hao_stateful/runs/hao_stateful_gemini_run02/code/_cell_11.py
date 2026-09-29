# Notice that the batches are E2L1-8 and L1-5! Let's check if the clusters are separated by batch/technology or biology!
# Let's check batch distribution across clusters or if batch correction (e.g. Harmony) is needed or if Leiden 0.5 without batch correction is driven by batch.
adata_proc.obs['batch'] = [x.split('_')[0] for x in adata_proc.obs_names]
adata_proc.obs['platform'] = ['Chromium' if x.startswith('E2') else 'Other' for x in adata_proc.obs['batch']]
print(pd.crosstab(adata_proc.obs['leiden_0.5'], adata_proc.obs['platform']))