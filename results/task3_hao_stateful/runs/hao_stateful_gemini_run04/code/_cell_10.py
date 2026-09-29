# Let's check batch distribution across clusters!
batch_crosstab = pd.crosstab(adata.obs['leiden_0.5'], adata.obs['batch'])
print(batch_crosstab)