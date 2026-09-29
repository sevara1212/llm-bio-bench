# Let's check expression of canonical markers across clusters:
sc.pl.dotplot(adata_filtered, marker_genes, groupby='leiden', show=False)
# Let's check mean expression of these markers per cluster
cluster_means = pd.DataFrame(index=groups, columns=marker_genes)
for g in groups:
    cells_in_g = adata_filtered.obs['leiden'] == g
    # adata_filtered.raw.X is log-normalized counts
    for m in marker_genes:
        if m in adata_filtered.raw.var_names:
            idx = adata_filtered.raw.var_names.get_loc(m)
            cluster_means.loc[g, m] = np.mean(adata_filtered.raw.X[cells_in_g, idx].toarray())
print(cluster_means.round(2))