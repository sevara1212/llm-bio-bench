# Let's inspect each cluster's top markers:
for group in sorted(adata.obs['leiden'].unique()):
    top_genes = [result['names'][group][i] for i in range(10)]
    print(f"Cluster {group}: {', '.join(top_genes)}")