for m in marker_genes:
    if m in adata_filtered.raw.var_names:
        idx = adata_filtered.raw.var_names.get_loc(m)
        vals = adata_filtered.raw.X[:, idx].toarray().flatten()
        for g in groups:
            cells_in_g = (adata_filtered.obs['leiden'] == g).values
            cluster_means.loc[g, m] = np.mean(vals[cells_in_g])
print(cluster_means.round(2))