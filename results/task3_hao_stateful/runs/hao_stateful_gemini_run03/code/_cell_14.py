# Let's inspect each cluster carefully:
# Cluster 0: RPL/RPS ribosomal genes, high CD3D (1.83), IL7R (2.17), CCR7 (0.85) -> Naive / Memory CD4+ T cell
# Cluster 1: MS4A1, CD79A, BANK1, RALGPS2 -> B cells
# Cluster 2: NKG7, CCL5, GZMH, CST7, CD8A (0.59), CD3D (1.98) -> Effector CD8+ T cells / cytotoxic T cells
# Cluster 3: IL7R, IL32, TPT1, LTB, CD3D (1.68), IL7R (2.99) -> CD4+ T cells (Memory CD4+ T cells)
# Cluster 4: STMN1, DEK, HMGB2, HMGN2, TUBA1B -> Proliferating cells / cycling lymphocytes
# Cluster 5: VCAN, CD14, BLVRB, FCN1, MNDA, S100A12, S100A8 -> CD14+ Monocytes
# Cluster 6: KLRB1, DUSP2, IL7R, GZMK, CD69, CD3D (1.84) -> CD8+ Memory / Effector Memory T cells (GZMK+ T cells)
# Cluster 7: LST1, AIF1, CFD, MS4A7, FTL, COTL1 -> FCGR3A+ / CD16+ Monocytes (non-classical monocytes)
# Cluster 8: KLRD1, GNLY, XCL2, CTSW, IL2RB, KLRF1, CD3D low (0.20) -> NK cells
# Cluster 9: HLA-DRA, HLA-DPA1, HLA-DPB1, CST3, CD14 low, HLA high -> Conventional Dendritic Cells (cDCs)
# Cluster 10: CCDC50, TCF4, PLD4, IRF7, LILRA4 -> Plasmacytoid Dendritic Cells (pDCs)
# Cluster 11: MZB1, HSP90B1, JCHAIN, TNFRSF17 (BCMA) -> Plasma cells / Plasmablasts
# Cluster 12: CDK6, SOX4, PRSS57 -> Progenitor cells / HSPCs
# Cluster 13: PF4, PPBP, CAVIN2, GNG11, TUBB1 -> Platelets / Megakaryocytes
# Cluster 14: PRF1, NKG7, CST7, CD247, KLRF1, KLRD1, GZMB, GNLY high (4.37), CD3D low (0.18) -> NK cells (CD56dim / cytotoxic NK)
# Cluster 15: GZMK, CD27, NUCB2, TCF7, CD3D (2.26) -> CD8+ T cells (Central memory / GZMK+)
# Cluster 16: CD74, HLA-DPA1, C1orf54, WDFY4, CPVL -> Dendritic cells (cDC1 / DC)
# Cluster 17: TNFRSF18 (GITR), TNFRSF4 (OX40), KLF2, CD3D/IL7R -> Regulatory T cells (Tregs) / Activated T cells
# Cluster 18: PPP1R14A, SAMHD1, PLD4, LILRA4 -> Plasmacytoid Dendritic Cells (pDCs)

# Let's check marker expression for cDC1 vs cDC2, Tregs (FOXP3), etc.
for g in ['CD8B', 'CD4', 'FCGR3A', 'MS4A7', 'FOXP3', 'CLEC9A', 'CLEC10A', 'CD1C', 'FCER1A']:
    if g in adata_proc.raw.var_names:
        print(f"{g}:", [round(float(adata_proc[adata_proc.obs['leiden_0.5'] == str(i), g].X.mean()), 2) for i in range(19)])