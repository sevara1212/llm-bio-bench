# Check marker genes for 0, 1, 4, 5, 6, 9 specifically
# Cluster 4 vs 5:
# Cluster 4: FCGR3A (CD16) high, LST1, CFD, SERPINA1 -> Non-classical / CD16+ Monocytes!
# Cluster 5: CD14, LYZ, CST3, CD1C -> Classical / CD14+ Monocytes!
# Cluster 6: CTSW, CD7, XCL1, XCL2, KLRD1, GNLY, TRDC -> NK cells (specifically CD56bright or NK/gdT cells)!
# Cluster 1: NKG7, CST7, GZMA, CCL5, PRF1, CD8A high, CD3D/E high -> Cytotoxic T cells / CD8+ T cells (or CD8+ T / NK)!
# Cluster 0: CD3D high, CD3E high, IL7R, CD4/CD8? Let's check CD4, CD8A, IL7R, CCR7, TCF7 in Cluster 0 vs 1 vs 12
check_genes = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'IL7R', 'CCR7', 'TCF7', 'LEF1', 'FOXP3', 'GZMK', 'NKG7', 'NCAM1', 'FCGR3A', 'CD14', 'MS4A7', 'CD34', 'SOX4', 'KIT', 'AVP']
for g in check_genes:
    if g in adata_proc.var_names:
        c0 = adata_proc[adata_proc.obs['leiden_0.3'] == '0', g].X.mean()
        c1 = adata_proc[adata_proc.obs['leiden_0.3'] == '1', g].X.mean()
        c6 = adata_proc[adata_proc.obs['leiden_0.3'] == '6', g].X.mean()
        c9 = adata_proc[adata_proc.obs['leiden_0.3'] == '9', g].X.mean()
        c12 = adata_proc[adata_proc.obs['leiden_0.3'] == '12', g].X.mean()
        print(f"{g:8s} Cl0: {c0:.2f} | Cl1: {c1:.2f} | Cl6: {c6:.2f} | Cl9: {c9:.2f} | Cl12: {c12:.2f}")