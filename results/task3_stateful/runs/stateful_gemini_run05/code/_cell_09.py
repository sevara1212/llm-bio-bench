sc.tl.rank_genes_groups(adata_filtered, 'leiden', method='t-test')
result = adata_filtered.uns['rank_genes_groups']
groups = result['names'].dtype.names
marker_df = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals_adj']}).head(5)
print(marker_df)