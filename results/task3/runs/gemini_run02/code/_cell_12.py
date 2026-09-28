# Let's inspect cluster sizes and marker expression across clusters
markers = ['IL7R', 'CD3D', 'CD4', 'CD8A', 'CD8B', 'S100A8', 'CD14', 'FCGR3A', 'MS4A7', 
           'MS4A1', 'CD79A', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP']

# Check mean expression per cluster
mean_df = pd.DataFrame(index=markers, columns=range(9))
for g in markers:
    if g in adata.raw.var_names:
        # extract expression from raw
        expr = adata.raw[:, g].X.toarray().flatten()
        for c in range(9):
            mean_df.loc[g, c] = expr[adata.obs['leiden'] == str(c)].mean()

print(mean_df)
print("\nCluster counts:")
print(adata.obs['leiden'].value_counts().sort_index())