# Let's inspect each cluster's top differentially expressed genes more thoroughly
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')
for cl in range(19):
    genes = [adata.uns['rank_genes_groups']['names'][str(cl)][i] for i in range(15)]
    scores = [adata.uns['rank_genes_groups']['scores'][str(cl)][i] for i in range(5)]
    print(f"Cluster {cl} (n={sum(adata.obs['leiden_0.5']==str(cl))}): {', '.join(genes)}")