# Let's map out each cluster:
# Cluster 0: Naive / Memory CD4+ T cell (CD3D+, IL7R+, CCR7+)
# Cluster 1: B cell (MS4A1+, CD79A+, BANK1+)
# Cluster 2: CD8+ T cell / Cytotoxic T cell (CD3D+, CD8A+, NKG7+, CCL5+, GZMH+)
# Cluster 3: CD4+ T cell (CD3D+, IL7R+, IL32+, LTB+)
# Cluster 4: Proliferating T/NK cell (STMN1+, HMGB2+, TUBA1B+)
# Cluster 5: CD14+ Monocyte (CD14+, VCAN+, FCN1+, S100A12+)
# Cluster 6: CD8+ T cell / Effector Memory T cell (CD3D+, IL7R+, GZMK+, KLRB1+)
# Cluster 7: FCGR3A+ Monocyte / CD16+ Monocyte (FCGR3A+, MS4A7+, LST1+, AIF1+)
# Cluster 8: NK cell (GNLY+, NKG7+, KLRD1+, KLRF1+, CD3D-)
# Cluster 9: Dendritic cell / cDC2 (FCER1A+, CD1C+, CLEC10A+, HLA-DRA+)
# Cluster 10: Plasmacytoid Dendritic cell / pDC (LILRA4+, CCDC50+, TCF4+, PLD4+, IRF7+)
# Cluster 11: Plasma cell (MZB1+, JCHAIN+, TNFRSF17+, HSP90B1+)
# Cluster 12: Hematopoietic stem and progenitor cell / HSPC (CDK6+, SOX4+, PRSS57+)
# Cluster 13: Platelet (PF4+, PPBP+, CAVIN2+, GNG11+)
# Cluster 14: NK cell (GNLY+, NKG7+, PRF1+, GZMB+, FCGR3A+, CD3D-)
# Cluster 15: CD8+ T cell (CD3D+, CD8B+, GZMK+, CD27+, TCF7+)
# Cluster 16: Dendritic cell / cDC1 (CLEC9A+, WDFY4+, C1orf54+, CD74+)
# Cluster 17: CD4+ T cell / Regulatory T cell (CD3D+, IL7R+, TNFRSF18+, TNFRSF4+)
# Cluster 18: Plasmacytoid Dendritic cell / pDC (LILRA4+, SAMHD1+, PLD4+, PPP1R14A+)

# Let's verify standard PBMC cell type naming:
cluster_to_celltype = {
    '0': 'CD4+ T cell',
    '1': 'B cell',
    '2': 'CD8+ T cell',
    '3': 'CD4+ T cell',
    '4': 'Proliferating T cell',
    '5': 'CD14+ Monocyte',
    '6': 'CD8+ T cell',
    '7': 'CD16+ Monocyte',
    '8': 'NK cell',
    '9': 'Dendritic cell',
    '10': 'Plasmacytoid dendritic cell',
    '11': 'Plasma cell',
    '12': 'Hematopoietic stem and progenitor cell',
    '13': 'Platelet',
    '14': 'NK cell',
    '15': 'CD8+ T cell',
    '16': 'Dendritic cell',
    '17': 'CD4+ T cell',
    '18': 'Plasmacytoid dendritic cell'
}

adata_proc.obs['cell_type'] = adata_proc.obs['leiden_0.5'].map(cluster_to_celltype)
print(adata_proc.obs['cell_type'].value_counts())

# Save labels.csv with columns barcode,cell_type
labels_df = pd.DataFrame({
    'barcode': adata_proc.obs_names,
    'cell_type': adata_proc.obs['cell_type'].values
})
labels_df.to_csv('labels.csv', index=False)
print("Saved labels.csv with shape:", labels_df.shape)
print("First 10 rows:")
print(labels_df.head(10))