# Let's check marker expression across clusters for res=1.0:
marker_genes = ['IL7R', 'CCR7', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP', 'CD3D']
for mg in marker_genes:
    if mg in adata.raw.var_names:
        print(f"Gene {mg} is present in adata.raw")
    else:
        print(f"Gene {mg} NOT present in adata.raw")

# Let's compute mean expression of these markers per cluster at res=1.0
cluster_means = pd.DataFrame(index=marker_genes, columns=range(9))
for cluster in range(9):
    cells = adata_processed.obs['leiden'] == str(cluster)
    for mg in marker_genes:
        idx = adata.raw.var_names.get_loc(mg)
        val = adata.raw.X[cells, idx].mean()
        cluster_means.loc[mg, cluster] = val
print(cluster_means.round(3))