# Let's inspect cluster 0, 1, 2, 3, 4, 5 more carefully:
# Let's check marker expression for well-known cell types:
# T cells: CD3D, CD3E, CD3G, CD4, CD8A, CD8B, IL7R
# Monocytes (CD14+ Monocytes, FCGR3A+ Monocytes): CD14, LYZ, CST3, MS4A7, FCGR3A
# NK cells: GNLY, NKG7, PRF1, NCAM1
# B cells: CD79A, CD79B, MS4A1 (CD20)
# Dendritic cells: FCER1A, CST3, HLA-DPA1/DPB1/DRA without B cell markers
# Platelets / Megakaryocytes: PPBP, PF4, GPX1
# What about cluster 0 vs 2?
# Let's check expression of CD3D, CD4, CD8A, GNLY, NKG7 across clusters:
marker_genes = ['CD3D', 'CD4', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'MS4A1', 'CD79A', 'CD14', 'FCGR3A', 'CST3', 'FCER1A', 'PPBP']
for g in marker_genes:
    if g in adata.raw.var_names:
        print(f"{g}:")
        mean_exp = [adata.raw[:, g].X[adata.obs['leiden'] == str(c)].mean() for c in range(6)]
        print("  ", [f"C{c}: {mean_exp[c]:.2f}" for c in range(6)])