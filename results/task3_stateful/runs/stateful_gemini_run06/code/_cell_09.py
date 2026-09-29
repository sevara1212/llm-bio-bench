genes_to_check = ['IL7R', 'CCR7', 'S100A4', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 'CD8B', 'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP', 'CD3D']
# compute mean expression per cluster in raw/log-normalized counts
mean_expr = pd.DataFrame(index=genes_to_check, columns=sorted(adata_scaled.obs['leiden'].unique()))
for g in genes_to_check:
    if g in adata.raw.var_names:
        for c in mean_expr.columns:
            cells = adata_scaled.obs['leiden'] == c
            mean_expr.loc[g, c] = adata.raw[cells, g].X.toarray().mean()

print(mean_expr.round(2))