# Let's inspect cluster 0, 11, 16 more closely
print("Cluster 0 top 20 genes:")
print(pd.DataFrame(adata.uns['rank_genes_groups']['names'])['0'].head(20).tolist())

print("\nCluster 11 top 20 genes:")
print(pd.DataFrame(adata.uns['rank_genes_groups']['names'])['11'].head(20).tolist())

print("\nCluster 16 top 20 genes:")
print(pd.DataFrame(adata.uns['rank_genes_groups']['names'])['16'].head(20).tolist())

# Also let's check CD4 vs CD8 vs NK vs B vs Mono in each cluster
genes = ['CD3D', 'CD4', 'CD8A', 'IL7R', 'CCR7', 'S100A4', 'FOXP3', 'NKG7', 'GNLY', 'NCAM1', 'FCGR3A', 'MS4A1', 'CD14', 'FCER1A', 'CST3', 'CLEC9A', 'LILRA4', 'TCF4', 'MZB1', 'PPBP', 'HBB', 'MKI67', 'CD34']
g_present = [g for g in genes if g in adata.raw.var_names]
sub_df = pd.DataFrame(index=[str(i) for i in range(19)])
for g in g_present:
    sub_df[g] = [adata.raw[adata.obs['leiden_0.5'] == str(i), g].X.mean() for i in range(19)]
print("\nExpression matrix:")
print(sub_df[['CD3D', 'CD4', 'CD8A', 'IL7R', 'NKG7', 'GNLY', 'MS4A1', 'CD14', 'FCGR3A', 'PPBP', 'HBB', 'MZB1', 'TCF4']])