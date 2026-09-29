# Notice min n_genes_by_counts is 499! Max pct_counts_mt is 15.06%!
# The dataset might already be pre-filtered or very high quality!
# Let's check where the dataset might have come from (e.g. 10k PBMC, or similar standard 10x dataset).
# Let's check barcode prefixes:
import scanpy as sc
adata = sc.read_h5ad('raw_counts.h5ad')
print(adata.obs_names[:20])
print(adata.obs_names.str.split('_').str[0].value_counts())