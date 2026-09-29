# Let's inspect T cell clusters (0, 1, 2, 8, 14, 16, 18, 4) and B cell clusters (3, 6):
# Cluster 3 vs 6:
# 3: IGKC, BANK1, MS4A1, CD79A
# 6: IGLC2, IGLC3, MS4A1, BANK1 -> Kappa vs Lambda B cells! Both are B cells (or Naive/Memory B cells).
# Cluster 1: FGFBP2, GZMH, GNLY, PRF1, CD8A (Effector CD8+ T cells / cytotoxic T cells)
# Cluster 2: GZMK, DUSP2, KLRB1, CD8A (Memory CD8+ T cells)
# Cluster 0: IL7R, TCF7, LEF1, CD3D (Naive CD4+ T cells)
# Cluster 16: GZMK, CD27, CD3D, TCF7, LCK (CD8+ / CD4+ central memory T cells)
# Cluster 8: KLRD1, GNLY, XCL1, XCL2, KLRF1, KLRC1 (NK cells - CD3D is negative!)
# Cluster 18: MKI67, CENPF, BIRC5, PRF1, NKG7 (Proliferating NK / T cells)
# Cluster 4: STMN1, DEK, PCNA, TYMS, CD8A, FCGR3A (Proliferating cytotoxic T / NK cells)
# Cluster 12: CD14+ Monocytes
# Cluster 7: FCGR3A+ (CD16+) Monocytes
# Cluster 5: Conventional Dendritic Cells (cDC2 / CD1c+ DC)
# Cluster 15: Conventional Dendritic Cells (cDC1 / CLEC9A+ DC)
# Cluster 9: Plasmacytoid Dendritic Cells (pDC)
# Cluster 10: Plasma cells
# Cluster 11: Hematopoietic Stem and Progenitor Cells (HSPC / CD34+ progenitors)
# Cluster 13: Platelets / Megakaryocytes
# Cluster 17: Erythrocytes / Red Blood Cells

# Let's verify cluster 8 CD3D:
print("Cluster 8 CD3D:", adata[adata.obs['leiden'] == '8', 'CD3D'].X.mean()) # 0.218 -> NK cells
print("Cluster 1 CD3D:", adata[adata.obs['leiden'] == '1', 'CD3D'].X.mean()) # 1.51 -> CD8+ Effector T cells
print("Cluster 0 CD3D:", adata[adata.obs['leiden'] == '0', 'CD3D'].X.mean()) # 1.81 -> CD4+ Naive T cells
print("Cluster 2 CD3D:", adata[adata.obs['leiden'] == '2', 'CD3D'].X.mean()) # 1.91 -> CD8+ Memory T cells
print("Cluster 14 CD3D:", adata[adata.obs['leiden'] == '14', 'CD3D'].X.mean()) # IL7R+, CD4+ Memory T cells

cluster_annotations = {
    '0': 'CD4+ Naive T cells',
    '1': 'CD8+ Effector T cells',
    '2': 'CD8+ Memory T cells',
    '3': 'B cells',
    '4': 'Proliferating T/NK cells',
    '5': 'Dendritic cells',
    '6': 'B cells',
    '7': 'CD16+ Monocytes',
    '8': 'NK cells',
    '9': 'Plasmacytoid dendritic cells',
    '10': 'Plasma cells',
    '11': 'HSPCs',
    '12': 'CD14+ Monocytes',
    '13': 'Platelets',
    '14': 'CD4+ Memory T cells',
    '15': 'Dendritic cells',
    '16': 'CD8+ T cells',
    '17': 'Erythrocytes',
    '18': 'Proliferating T/NK cells'
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)
df_out = pd.DataFrame({'barcode': adata.obs_names, 'cell_type': adata.obs['cell_type']})
df_out.to_csv('labels.csv', index=False)
print("Saved labels.csv with shape:", df_out.shape)
print(df_out['cell_type'].value_counts())
print(df_out.head())