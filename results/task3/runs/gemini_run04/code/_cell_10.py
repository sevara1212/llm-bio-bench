# Let's inspect res=1.0
sc.tl.rank_genes_groups(adata_hvg, 'leiden_1.0', method='t-test')
result = adata_hvg.uns['rank_genes_groups']
print("--- Resolution 1.0 ---")
for g in result['names'].dtype.names:
    top_genes = list(result['names'][g][:8])
    print(f"Cluster {g} (n={sum(adata_hvg.obs['leiden_1.0']==g)}): {top_genes}")