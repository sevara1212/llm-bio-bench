import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('processed.h5ad')

# Print details of clusters 8, 14, 5, 7, 9, 10, 11, 12, 16, 18
cl_markers = {
    '0': 'CD4+ Naive T cells', # CD3D+, IL7R+, CCR7 high, RPL/RPS
    '1': 'B cells',            # MS4A1+, CD79A+, BANK1+
    '2': 'CD8+ Effector T cells', # CD3D+, CD8A+, GZMH+, CCL5+, NKG7+
    '3': 'CD4+ Memory T cells', # CD3D+, CD4+, IL7R+, LTB+, S100A4+
    '4': 'Proliferating T cells', # STMN1+, HMGB2+, TUBA1B+, CD3D+
    '5': 'CD14+ Monocytes',    # CD14+, VCAN+, FCN1+, S100A8/9/12+
    '6': 'MAIT cells',         # KLRB1+, SLC4A10+, RORA+, DUSP2+
    '7': 'FCGR3A+ Monocytes',  # FCGR3A/CD16+, LST1+, MS4A7+, AIF1+
    '8': 'NK cells',           # GNLY+, NKG7+, KLRD1+, KLRF1+, CD3D-
    '9': 'Conventional Dendritic Cells', # HLA-DRA/DP/DQ+, CST3+, LYZ+, CD1C+
    '10': 'Plasmacytoid Dendritic Cells', # TCF4+, CCDC50+, LILRA4+, PLD4+, IRF7/8+
    '11': 'Plasma cells',      # MZB1+, JCHAIN+, SDC1/TNFRSF17+
    '12': 'Hematopoietic stem and progenitor cells', # CDK6+, SOX4+, PRSS57+
    '13': 'Platelets',         # PPBP+, PF4+, NRGN+, TUBB1+
    '14': 'NK cells',          # FCGR3A+, PRF1+, NKG7+, KLRD1+, CD3D-
    '15': 'CD8+ Central Memory T cells', # CD3D+, CD8B+, GZMK+, TCF7+, CD27+
    '16': 'Conventional Dendritic Cells', # CLEC9A+, WDFY4+, CD74+, HLA-DP/DQ/DR+
    '17': 'Gamma-delta T cells', # TRDC+, TRGC1+, TNFRSF18+, KLRB1+
    '18': 'Plasmacytoid Dendritic Cells', # LILRA4+, PLD4+, SAMHD1+, PPP1R14A+, CST3+
}

# Standard canonical PBMC mapping:
# 0: CD4+ Naive T cells (or CD4+ T cells)
# 1: B cells
# 2: CD8+ Effector T cells (or CD8+ T cells)
# 3: CD4+ Memory T cells (or CD4+ T cells)
# 4: Proliferating T cells (or Proliferating cells)
# 5: CD14+ Monocytes (or Monocytes)
# 6: MAIT cells (or CD8+ T cells / MAIT)
# 7: CD16+ Monocytes (or Monocytes)
# 8: NK cells
# 9: Dendritic cells (or Conventional dendritic cells)
# 10: Plasmacytoid dendritic cells
# 11: Plasma cells
# 12: Progenitor cells (or HSPC)
# 13: Platelets
# 14: NK cells (or CD16+ NK cells)
# 15: CD8+ Memory T cells (or CD8+ T cells)
# 16: Conventional dendritic cells (or Dendritic cells)
# 17: Gamma-delta T cells
# 18: Plasmacytoid dendritic cells

# Let's map directly:
cell_type_mapping = {
    '0': 'CD4+ Naive T cells',
    '1': 'B cells',
    '2': 'CD8+ Effector T cells',
    '3': 'CD4+ Memory T cells',
    '4': 'Proliferating T cells',
    '5': 'CD14+ Monocytes',
    '6': 'MAIT cells',
    '7': 'FCGR3A+ Monocytes',
    '8': 'NK cells',
    '9': 'Conventional Dendritic Cells',
    '10': 'Plasmacytoid Dendritic Cells',
    '11': 'Plasma cells',
    '12': 'Hematopoietic Stem and Progenitor Cells',
    '13': 'Platelets',
    '14': 'NK cells',
    '15': 'CD8+ Memory T cells',
    '16': 'Conventional Dendritic Cells',
    '17': 'Gamma-delta T cells',
    '18': 'Plasmacytoid Dendritic Cells'
}

adata.obs['cell_type'] = adata.obs['leiden_0.5'].map(cell_type_mapping)

# Ensure all 7841 raw barcodes are mapped
raw_adata = sc.read_h5ad('raw_counts.h5ad')
labels_df = pd.DataFrame(index=raw_adata.obs_names)
labels_df['barcode'] = raw_adata.obs_names
labels_df['cell_type'] = adata.obs.loc[raw_adata.obs_names, 'cell_type']

labels_df.to_csv('labels.csv', index=False)
print(labels_df['cell_type'].value_counts())
print("Null count:", labels_df['cell_type'].isnull().sum())
print("Shape:", labels_df.shape)
print("First 10 rows:")
print(labels_df.head(10))