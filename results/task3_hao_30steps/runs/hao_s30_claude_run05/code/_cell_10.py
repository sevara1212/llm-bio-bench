import scanpy as sc

adata = sc.read_h5ad('full_with_clusters.h5ad')

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# Print top markers for each cluster
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
for g in groups:
    top_genes = [result['names'][g][i] for i in range(10)]
    print(f"Cluster {g}: {top_genes}")