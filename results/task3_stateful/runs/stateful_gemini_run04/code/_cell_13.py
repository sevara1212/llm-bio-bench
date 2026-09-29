# Let's inspect res=0.8 clusters and mean expressions of key PBMC marker genes:
sc.tl.leiden(adata_proc, resolution=0.8, random_state=0)
adata_filtered.obs['leiden'] = adata_proc.obs['leiden']

marker_genes = ['IL7R', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP', 'CD3D']
df_expr = pd.DataFrame(index=sorted(adata_filtered.obs['leiden'].unique()))
for g in marker_genes:
    # get expression from raw
    g_idx = adata_filtered.raw.var_names.get_loc(g)
    expr = adata_filtered.raw.X[:, g_idx].toarray().flatten()
    df_expr[g] = [expr[adata_filtered.obs['leiden'] == cl].mean() for cl in df_expr.index]

print(df_expr.round(2))