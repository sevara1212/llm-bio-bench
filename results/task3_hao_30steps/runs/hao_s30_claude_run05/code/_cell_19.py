import scanpy as sc

adata = sc.read_h5ad('full_with_clusters.h5ad')
result = adata.uns.get('rank_genes_groups')
# recompute if needed, but let's check obs leiden categories match earlier count
print(adata.obs['leiden'].value_counts())