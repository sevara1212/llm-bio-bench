import pandas as pd
import scanpy as sc

adata = sc.read_h5ad('processed.h5ad')

# Cell type mappings based on Leiden 0.5 clusters:
# 0: Naive CD4+ T cell (LEF1, TCF7, CCR7, CD3D, IL7R)
# 1: Memory CD4+ T cell (IL7R, LTB, S100A4, CD3D)
# 2: CD8+ Effector T cell (CCL5, NKG7, GZMH, CST7, CD8A/B, CD3D)
# 3: Memory B cell (MS4A1, BANK1, CD79A, CD37, CD27)
# 4: Proliferating / Cycling T cell (STMN1, HMGB2, H2AZ1, TUBA1B, CD3D)
# 5: MAIT / Effector Memory T cell (KLRB1, DUSP2, GZMK, IL7R, CD3D)
# 6: CD14+ Monocyte (CD14, VCAN, S100A9, S100A8, LYZ)
# 7: NK cell (PRF1, NKG7, GZMB, KLRD1, GNLY, NCAM1)
# 8: CD16+ Monocyte (FCGR3A/CD16, LST1, AIF1, MS4A7, CFD)
# 9: NK cell (CD56bright / Adaptive NK: KLRD1, GNLY, XCL1, XCL2, KLRF1, KLRC1)
# 10: Naive B cell (TCL1A, IGHD, IGHM, CD79A, MS4A1)
# 11: Plasmacytoid Dendritic Cell (pDC) (LILRA4, IL3RA, TCF4, CCDC50, PLD4)
# 12: Hematopoietic Stem / Progenitor Cell (HSPC) (CD34, CDK6, SOX4, PRSS57)
# 13: Plasma cell / Plasmablast (MZB1, JCHAIN, SDC1, TNFRSF17, SEC11C)
# 14: Platelet / Megakaryocyte (PPBP, PF4, CAVIN2, GNG11)
# 15: Conventional Dendritic Cell 2 (cDC2) (FCER1A, HLA-DRA, CST3, CD1C)
# 16: CD8+ Central Memory T cell (GZMK, CD27, TCF7, CD8A, CD3D)
# 17: Conventional Dendritic Cell 1 (cDC1) (CLEC9A, WDFY4, CPVL, C1orf54, CST3)
# 18: Regulatory T cell / Activated T cell (TNFRSF18, TNFRSF4, KLRB1, IL7R, CD3D)
# 19: Erythroid cell (HBB, HBA1, HBA2, ALAS2, SLC25A37)
# 20: Plasmacytoid Dendritic Cell (pDC) (LILRA4, SAMHD1, PLD4, PPP1R14A)

cluster2celltype = {
    '0': 'Naive CD4+ T cell',
    '1': 'Memory CD4+ T cell',
    '2': 'CD8+ Effector T cell',
    '3': 'Memory B cell',
    '4': 'Proliferating T cell',
    '5': 'MAIT / Effector Memory T cell',
    '6': 'CD14+ Monocyte',
    '7': 'NK cell',
    '8': 'CD16+ Monocyte',
    '9': 'NK cell',
    '10': 'Naive B cell',
    '11': 'Plasmacytoid Dendritic Cell',
    '12': 'HSPC',
    '13': 'Plasma cell',
    '14': 'Platelet',
    '15': 'cDC2',
    '16': 'CD8+ Central Memory T cell',
    '17': 'cDC1',
    '18': 'Treg / Activated T cell',
    '19': 'Erythroid cell',
    '20': 'Plasmacytoid Dendritic Cell'
}

raw_adata = sc.read_h5ad('raw_counts.h5ad')
df_labels = pd.DataFrame({
    'barcode': raw_adata.obs_names,
    'cell_type': raw_adata.obs_names.map(adata.obs['leiden_0.5'].map(cluster2celltype))
})

print(df_labels.head())
print("Null count:", df_labels.isnull().sum().to_dict())
print("Cell type distribution:")
print(df_labels['cell_type'].value_counts())

df_labels.to_csv('labels.csv', index=False)
print("Saved labels.csv successfully!")