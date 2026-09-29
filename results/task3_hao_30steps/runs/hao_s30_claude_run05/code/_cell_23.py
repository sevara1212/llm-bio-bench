import scanpy as sc
adata = sc.read_h5ad('full_with_clusters.h5ad')

result = adata.uns['rank_genes_groups']
for g in ['18','20']:
    top_genes = [result['names'][g][i] for i in range(15)]
    print(f"Cluster {g}: {top_genes}")