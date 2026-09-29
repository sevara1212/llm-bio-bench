# In Scanpy tutorial:
# leiden with default resolution (which was 1.0 or 0.8 in older versions? In tutorial: sc.tl.leiden(adata) which defaults to resolution=1.0)
# Let's test resolution=0.8, 1.0
sc.tl.leiden(adata, resolution=0.8)
print("Resolution 0.8 clusters:", adata.obs['leiden'].value_counts())

# Let's find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
marker_df = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals_adj']}).head(10)
print(marker_df)