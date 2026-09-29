# Top DE genes for cluster 0, 1, 12, 6, 9
sc.tl.rank_genes_groups(adata_proc, groupby='leiden_0.3', method='wilcoxon')
for cl in ['0', '1', '4', '5', '6', '9', '12']:
    genes = adata_proc.uns['rank_genes_groups']['names'][cl][:10]
    scores = adata_proc.uns['rank_genes_groups']['scores'][cl][:10]
    print(f"Cluster {cl}: {list(zip(genes, np.round(scores, 2)))}")