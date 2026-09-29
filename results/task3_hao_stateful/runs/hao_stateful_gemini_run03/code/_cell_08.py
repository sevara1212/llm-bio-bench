# Let's check barcode prefixes
prefixes = adata.obs_names.map(lambda x: x.split('_')[0]).value_counts()
print(prefixes)