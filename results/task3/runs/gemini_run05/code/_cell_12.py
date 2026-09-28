# Let's inspect all clusters in leiden_0.8 and leiden_1.0
sc.tl.rank_genes_groups(adata, groupby='leiden_0.8', method='wilcoxon')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
for g in groups:
    print(f"Cluster {g} (N={sum(adata.obs['leiden_0.8'] == g)}):", [result['names'][g][i] for i in range(12)])

sc.tl.rank_genes_groups(adata, groupby='leiden_1.0', method='wilcoxon')
result1 = adata.uns['rank_genes_groups']
groups1 = result1['names'].dtype.names
print("\n--- Res 1.0 ---")
for g in groups1:
    print(f"Cluster {g} (N={sum(adata.obs['leiden_1.0'] == g)}):", [result1['names'][g][i] for i in range(12)])