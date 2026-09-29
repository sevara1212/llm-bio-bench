# Let's inspect cluster 0, 12, 17, 15, 4 etc. more carefully
# What are the major PBMC cell types?
# T cells: CD4+ T cells, CD8+ T cells, NK cells, B cells, Plasma cells, Monocytes (CD14+ Monocytes, FCGR3A+/CD16+ Monocytes),
# Dendritic cells (cDC1, cDC2 / myeloid DC, pDC), Megakaryocytes/platelets, Proliferating cells / Hematopoietic progenitors.

# Let's check marker expression for clusters:
# Cluster 0: RPL/RPS ribosomal genes high, CD3D=1.83, CD8B=0.50, IL7R?
# Let's see what T cell markers are in cluster 0, 3, 2, 6, 15, 17
t_markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'IL7R', 'CCR7', 'S100A4', 'FOXP3', 'GZMA', 'GZMB', 'GZMK', 'NKG7', 'NCAM1', 'FCGR3A', 'TRAC', 'TRDC']
t_present = [g for g in t_markers if g in adata.raw.var_names]

df_t = pd.DataFrame(index=[0, 2, 3, 4, 6, 8, 12, 14, 15, 17])
for g in t_present:
    g_idx = adata.raw.var_names.get_loc(g)
    expr = adata.raw.X[:, g_idx].toarray().flatten()
    df_t[g] = [expr[adata.obs['leiden_0.5'] == str(cl)].mean() for cl in df_t.index]

print(df_t.round(2))