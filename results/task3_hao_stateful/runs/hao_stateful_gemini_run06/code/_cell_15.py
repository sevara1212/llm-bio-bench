# Let's check cluster 11 markers and what kind of cells they are
print("Cluster 11 more markers:")
print(list(markers['11'])[:15])

# Let's map clusters to well-recognized PBMC cell types:
# 0: Naive / Memory CD4+ T cells (IL7R, TCF7, LEF1, CD3D, LTB)
# 1: Effector CD8+ T cells / cytotoxic T cells (CD3D, CD8A, GZMB, GZMH, PRF1, FGFBP2)
# 2: Memory / Effector Memory CD8+ T cells (GZMK+, KLRB1+, CD3D+, IL7R+)
# 3: B cells (MS4A1, CD79A, CD79B, BANK1, IGKC)
# 4: Proliferating T cells (STMN1, PCNA, TYMS, MKI67, CD3D)
# 5: Conventional Dendritic Cells / cDC2 (HLA-DRA, HLA-DPB1, CD1C, CST3)
# 6: B cells (MS4A1, CD79A, CD79B, IGLC2)
# 7: Non-classical monocytes (FCGR3A / CD16+, LST1, MS4A7, AIF1)
# 8: NK cells (KLRD1, KLRF1, GNLY, NCAM1, CTSW, CD3D negative)
# 9: Plasmacytoid Dendritic Cells / pDC (TCF4, IL3RA, CCDC50, PLD4)
# 10: Plasma cells (MZB1, JCHAIN, TNFRSF17, SDC1)
# 11: Hematopoietic stem/progenitor cells (HSPCs) / Progenitor cells (CDK6, SOX4, PRSS57)
# 12: Classical monocytes (CD14+, S100A8, S100A9, VCAN, FCN1)
# 13: Platelets / Megakaryocytes (PPBP, PF4, CAVIN2, NRGN, TUBB1)
# 14: Activated CD4+ T cells / T helper cells (IL7R, CD69, TNFRSF4, TNFRSF18, KLRB1)
# 15: Conventional Dendritic Cells (cDC1) (CLEC9A, WDFY4, SNX3)
# 16: Central Memory / Effector Memory CD8+ T cells (GZMK, CD8B, CD27, CD3D)
# 17: Erythrocytes / Erythroid cells (HBB, HBA1, HBA2, SLC25A37)
# 18: Proliferating NK/T cells (MKI67, CENPF, BIRC5, GNLY, NKG7)

cluster_names = {
    '0': 'CD4+ T cells',
    '1': 'Cytotoxic CD8+ T cells',
    '2': 'Effector Memory CD8+ T cells',
    '3': 'B cells',
    '4': 'Proliferating T cells',
    '5': 'cDC2',
    '6': 'B cells',
    '7': 'Non-classical monocytes',
    '8': 'NK cells',
    '9': 'Plasmacytoid dendritic cells',
    '10': 'Plasma cells',
    '11': 'HSPCs',
    '12': 'Classical monocytes',
    '13': 'Platelets',
    '14': 'Activated CD4+ T cells',
    '15': 'cDC1',
    '16': 'Memory CD8+ T cells',
    '17': 'Erythrocytes',
    '18': 'Proliferating NK cells'
}

adata_proc.obs['cell_type'] = adata_proc.obs['leiden_0.5'].map(cluster_names)
out_df = pd.DataFrame({
    'barcode': adata_proc.obs_names,
    'cell_type': adata_proc.obs['cell_type']
})

out_df.to_csv('labels.csv', index=False)
print("labels.csv saved. Head:")
print(out_df.head(10))
print("Value counts:")
print(out_df['cell_type'].value_counts())
print("Total rows:", len(out_df))