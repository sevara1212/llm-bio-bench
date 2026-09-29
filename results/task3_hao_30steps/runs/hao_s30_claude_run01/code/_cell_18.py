import scanpy as sc
adata = sc.read_h5ad('full_with_clusters.h5ad')

# look at top DE genes again for cluster 20
sc.tl.rank_genes_groups(adata, 'leiden', groups=['20'], method='wilcoxon')
res = adata.uns['rank_genes_groups']
print(list(res['names']['20'][:30]))