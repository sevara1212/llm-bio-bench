# Let's inspect cluster 0 more closely: why did RPL30 appear if non_ribo_mt filtered var_names? Ah, adata_sub was sliced from adata_proc which already had rank_genes_groups or maybe case-sensitivity.
# Let's look at Cluster 0 in adata_proc:
# Mean expression of CD3D is 1.83, CD3E 1.60, IL7R 2.17. That's CD4 T cells!
# Let's write a comprehensive script to verify markers for all clusters, assign labels, check and save labels.csv!

# Cluster 0: CD3D+, IL7R+ -> CD4 T cell
# Cluster 1: MS4A1, CD79A, BANK1 -> B cell
# Cluster 2: NKG7, CCL5, GZMH, CD8A/B -> CD8 T cell
# Cluster 3: IL7R, IL32, LTB, CD3D -> CD4 T cell / Naive CD4 T cell
# Cluster 4: STMN1, DEK, HMGB2, proliferating / dividing cells
# Cluster 5: VCAN, CD14, FCN1 -> CD14+ Monocyte
# Cluster 6: KLRB1, IL7R, GZMK, CD3D -> CD4/CD8 Memory T cell / T cell
# Cluster 7: LST1, MS4A7, FCGR3A -> FCGR3A+ / Non-classical Monocyte
# Cluster 8: GNLY, KLRD1, XCL2, NKG7 -> NK cell
# Cluster 9: HLA-DRA, HLA-DPA1, FCER1A, CST3 -> Dendritic cell (cDC)
# Cluster 10: TCF4, PLD4, LILRA4 -> Plasmacytoid Dendritic cell (pDC)
# Cluster 11: MZB1, JCHAIN -> Plasma cell
# Cluster 12: B cell progenitor / pre-B or B cell (CD79A+, CD19+, CDK6+, SOX4+)
# Cluster 13: PF4, PPBP -> Megakaryocyte / Platelet
# Cluster 14: PRF1, NKG7, CD247, KLRF1 -> NK cell / Cytotoxic T/NK
# Cluster 15: CD3D, CD27, TCF7, GZMK -> CD8 T cell / Central Memory T cell
# Cluster 16: CD74, HLA-DRA, CLEC9A -> Dendritic cell
# Cluster 17: TNFRSF18 (GITR), TNFRSF4 (OX40), Treg / Activated T cell
# Cluster 18: LILRA4, PLD4, PPP1R14A -> Plasmacytoid Dendritic cell (pDC)

# Let's check marker expression across major PBMC cell types:
# T cells (CD4+ T cells, CD8+ T cells, NK cells), B cells, Monocytes (CD14+, FCGR3A+), Dendritic cells, Platelets.
# Let's map clusters cleanly and save labels.csv.

cluster_map = {
    0: 'CD4 T cells',
    1: 'B cells',
    2: 'CD8 T cells',
    3: 'CD4 T cells',
    4: 'Proliferating T cells',
    5: 'CD14+ Monocytes',
    6: 'CD4 T cells',
    7: 'FCGR3A+ Monocytes',
    8: 'NK cells',
    9: 'Dendritic cells',
    10: 'Plasmacytoid Dendritic cells',
    11: 'Plasma cells',
    12: 'B cells',
    13: 'Platelets',
    14: 'NK cells',
    15: 'CD8 T cells',
    16: 'Dendritic cells',
    17: 'CD4 T cells',
    18: 'Plasmacytoid Dendritic cells'
}

adata_proc.obs['cell_type'] = adata_proc.obs['leiden_0.5'].astype(int).map(cluster_map)
df_out = pd.DataFrame({'barcode': adata_proc.obs_names, 'cell_type': adata_proc.obs['cell_type']})
df_out.to_csv('labels.csv', index=False)
print("labels.csv written successfully! Head:")
print(df_out.head(10))
print("\nCell type value counts:")
print(df_out['cell_type'].value_counts())
print("\nNull count:", df_out['cell_type'].isnull().sum())
print("Total rows:", len(df_out))