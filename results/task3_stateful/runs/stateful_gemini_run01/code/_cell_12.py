cells_arr = (adata_processed.obs['leiden'] == '0').values
cluster_means = pd.DataFrame(index=marker_genes, columns=[str(i) for i in range(9)])
for cluster in range(9):
    cells = (adata_processed.obs['leiden'] == str(cluster)).values
    for mg in marker_genes:
        idx = adata.raw.var_names.get_loc(mg)
        val = np.mean(adata.raw.X[cells, idx].toarray())
        cluster_means.loc[mg, str(cluster)] = val
print(cluster_means.round(3))