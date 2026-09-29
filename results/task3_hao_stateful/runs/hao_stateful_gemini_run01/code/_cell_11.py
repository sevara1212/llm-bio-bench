# Let's inspect each cluster more closely with canonical PBMC markers:
# CD3D, CD3E, CD4, CD8A, CD8B, NCAM1 (CD56), NKG7, GNLY, MS4A1 (CD20), CD19, CD14, FCGR3A (CD16),
# LYZ, HLA-DRA, CD1C, CLEC9A, CLEC4C/IL3RA/LILRA4 (pDC), PPBP/PF4 (Megakaryocytes/platelets), MZB1/TNFRSF17 (plasma cells)
canonical = ['CD3D', 'CD4', 'CD8A', 'CD8B', 'NCAM1', 'NKG7', 'GNLY', 'MS4A1', 'CD19', 'CD14', 'FCGR3A', 'LYZ', 'HLA-DRA', 'CD1C', 'CLEC9A', 'LILRA4', 'PPBP', 'MZB1', 'FOXP3']
canonical_present = [g for g in canonical if g in adata.raw.var_names]
sc.pl.dotplot(adata, canonical_present, groupby='leiden_0.5', save='_markers_0.5.png')

# Let's also check mean expression across clusters for these canonical markers
mean_exp = pd.DataFrame(index=adata.obs['leiden_0.5'].cat.categories)
for g in canonical_present:
    # get expression from raw
    g_idx = adata.raw.var_names.get_loc(g)
    expr = adata.raw.X[:, g_idx]
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    mean_exp[g] = [expr[adata.obs['leiden_0.5'] == cl].mean() for cl in adata.obs['leiden_0.5'].cat.categories]

print(mean_exp.round(2))