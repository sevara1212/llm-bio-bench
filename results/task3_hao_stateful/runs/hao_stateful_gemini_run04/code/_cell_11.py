# Let's run rank_genes_groups on leiden_0.5 and leiden_0.8
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')
marker_df_05 = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
print("Top 10 markers per cluster (leiden 0.5):")
print(marker_df_05)