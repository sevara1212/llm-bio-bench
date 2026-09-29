# Check cluster 0 and cluster 6 markers
sc.tl.rank_genes_groups(adata_processed, 'leiden', groups=['0', '6'], reference='rest', method='wilcoxon')
print("Cluster 0 top genes:", adata_processed.uns['rank_genes_groups']['names']['0'][:10])
print("Cluster 6 top genes:", adata_processed.uns['rank_genes_groups']['names']['6'][:10])

# Also compare 0 vs 6 directly
sc.tl.rank_genes_groups(adata_processed, 'leiden', groups=['6'], reference='0', method='wilcoxon')
print("Cluster 6 vs 0 top genes:", adata_processed.uns['rank_genes_groups']['names']['6'][:10])