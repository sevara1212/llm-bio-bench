# Let's inspect all clusters and their top markers in detail
result = adata.uns['rank_genes_groups']
for group in sorted(groups, key=lambda x: int(x)):
    print(f"\n--- Cluster {group} ---")
    top_genes = result['names'][group][:15]
    top_scores = result['scores'][group][:15]
    print(list(zip(top_genes, [round(s, 2) for s in top_scores])))