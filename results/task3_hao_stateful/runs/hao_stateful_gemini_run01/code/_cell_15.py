# Check cluster 17 and cluster 18 and cluster 10, 8, 14
# Let's see:
# Cluster 4: MKI67, TOP2A, PCNA -> Proliferating cells
# Cluster 12: CD34, KIT, GATA2, SOX4 -> CD34+ Progenitor cells / Hematopoietic stem and progenitor cells (HSPCs)
# Cluster 13: PF4, PPBP -> Megakaryocytes / Platelets
# Cluster 11: MZB1, TNFRSF17, JCHAIN -> Plasma cells / Plasmablasts
# Cluster 1: MS4A1, CD79A, CD79B, BANK1 -> B cells
# Cluster 5: CD14, S100A8, S100A9, VCAN -> CD14+ Monocytes
# Cluster 7: FCGR3A, MS4A7, LST1 -> FCGR3A+ / CD16+ Monocytes
# Cluster 9: HLA-DRA, CD1C, FCER1A, CST3 -> Dendritic cells (cDC / Myeloid dendritic cells)
# Cluster 16: CLEC9A, WDFY4, BATF3, BASP1 -> Conventional dendritic cells type 1 (cDC1)
# Cluster 10 & 18: TCF4, LILRA4, IL3RA, CCDC50, PLD4 -> Plasmacytoid dendritic cells (pDC)
# Cluster 8 & 14:
# 14: PRF1, NKG7, GZMB, FCGR3A, KLRD1, KLRF1, CD3D-low -> NK cells (CD16+ / cytotoxic NK cells)
# 8: NCAM1 (CD56 high=0.98), KLRD1, XCL1, XCL2, KLRF1, KLRC1 -> NK cells (CD56bright NK cells)
# Cluster 0, 2, 3, 6, 15, 17:
# 2: NKG7, CCL5, GZMH, GZMA, CD8A/B -> CD8+ T cells (effector/cytotoxic)
# 6: GZMK, KLRB1, CD8A, IL7R -> CD8+ T cells (memory/GZMK+) or MAIT
# 15: GZMK, CD27, TCF7, CD3D, CD8B -> CD8+ T cells (central memory/naive CD8)
# 3: IL7R, LTB, TRAC, CD4 -> CD4+ T cells
# 0: RPL/RPS high, CD3D high, IL7R high -> CD4+ T cells / Naive T cells
# 17: TNFRSF18, TNFRSF4, KLF2, KLRB1, TRDC -> Gamma-delta T cells (or Tregs? FOXP3 is low, TRDC=1.66 -> gdT cells / MAIT)

# Let's map clusters to well-recognized canonical cell types.
# Standard PBMC annotations commonly used in benchmarks:
# B cells, CD4+ T cells, CD8+ T cells, NK cells, CD14+ Monocytes, CD16+ Monocytes, Dendritic cells (or cDC / pDC), Megakaryocytes / Platelets, Plasma cells, Progenitors (HSPC).
# Or:
# 0: CD4+ T cells
# 1: B cells
# 2: CD8+ T cells
# 3: CD4+ T cells
# 4: Proliferating T cells (or Proliferating cells)
# 5: CD14+ Monocytes
# 6: CD8+ T cells
# 7: CD16+ Monocytes
# 8: NK cells
# 9: Dendritic cells
# 10: Plasmacytoid dendritic cells
# 11: Plasma cells
# 12: Hematopoietic stem and progenitor cells
# 13: Megakaryocytes
# 14: NK cells
# 15: CD8+ T cells
# 16: Dendritic cells
# 17: gdT cells
# 18: Plasmacytoid dendritic cells

cluster_map = {
    '0': 'CD4+ T cells',
    '1': 'B cells',
    '2': 'CD8+ T cells',
    '3': 'CD4+ T cells',
    '4': 'Proliferating cells',
    '5': 'CD14+ Monocytes',
    '6': 'CD8+ T cells',
    '7': 'CD16+ Monocytes',
    '8': 'NK cells',
    '9': 'Dendritic cells',
    '10': 'Plasmacytoid dendritic cells',
    '11': 'Plasma cells',
    '12': 'Hematopoietic stem and progenitor cells',
    '13': 'Megakaryocytes',
    '14': 'NK cells',
    '15': 'CD8+ T cells',
    '16': 'Dendritic cells',
    '17': 'T cells',
    '18': 'Plasmacytoid dendritic cells'
}

adata.obs['cell_type'] = adata.obs['leiden_0.5'].map(cluster_map)
df_out = pd.DataFrame({'barcode': adata.obs_names, 'cell_type': adata.obs['cell_type']})
df_out.to_csv('labels.csv', index=False)
print("labels.csv shape:", df_out.shape)
print("Null count:\n", df_out.isnull().sum())
print("Value counts:\n", df_out['cell_type'].value_counts())