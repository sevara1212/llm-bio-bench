# Let's check barcode names to see if there's batch/sample info
sample_prefixes = adata.obs_names.map(lambda x: x.split('_')[0]).value_counts()
print(sample_prefixes)