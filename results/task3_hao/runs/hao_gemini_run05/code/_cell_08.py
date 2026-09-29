import scanpy as sc

adata = sc.read_h5ad('raw_counts.h5ad')
print("Barcodes head:", adata.obs_names[:5].tolist())
print("Total cells:", len(adata.obs_names))

# Let's see if there are any batch/sample indicators in barcode prefix
samples = [b.split('_')[0] for b in adata.obs_names]
import collections
print("Sample prefixes:", collections.Counter(samples))