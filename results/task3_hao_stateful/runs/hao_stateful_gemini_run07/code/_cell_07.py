# Let's inspect the cells in raw_counts.h5ad
# Could this be from a public dataset or benchmark, like PBMC 10k or PBMC 3k or 8k, or multiome?
# Look at barcode names: "L3_ACCGTTCCACATATGC, E2L7_GACCAATCAGAGCTAG..."
# Let's check prefixes in barcodes:
barcodes = adata.obs_names
prefixes = [b.split('_')[0] for b in barcodes if '_' in b]
print("Prefix counts:", pd.Series(prefixes).value_counts())