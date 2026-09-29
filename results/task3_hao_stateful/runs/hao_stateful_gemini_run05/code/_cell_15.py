# Define cluster mapping
cluster_map = {
    '0': 'CD4+ T cells',            # High IL7R, CCR7, TCF7, LEF1, CD3D/E
    '1': 'CD8+ T cells',            # High CD8A, NKG7, CCL5, CST7, GZMA, PRF1, CD3D/E
    '2': 'B cells',                 # MS4A1 (CD20), CD79A, CD79B, BANK1, HLA-DRA
    '3': 'Proliferating cells',     # MKI67, TOP2A, STMN1, TUBA1B, HMGB2
    '4': 'CD16+ Monocytes',         # FCGR3A (CD16), MS4A7, LST1, FTL, AIF1
    '5': 'CD14+ Monocytes',         # CD14, LYZ, CST3, HLA-DRA
    '6': 'NK cells',                # NCAM1 (CD56), KLRD1, GNLY, XCL1, XCL2, CTSW
    '7': 'Plasmacytoid Dendritic cells', # TCF4, LILRA4, CLEC4C, IL3RA, CCDC50, IRF7
    '8': 'Plasma cells',            # MZB1, SDC1, JCHAIN, TNFRSF17, SEC11C
    '9': 'HSPC',                    # CD34, SOX4, AVP, CDK6 (hematopoietic stem/progenitor cells)
    '10': 'Platelets',              # PPBP, PF4, CAVIN2, NRGN, TUBB1, GNG11
    '11': 'Conventional Dendritic cells', # CLEC9A, WDFY4, CPVL, HLA-DRA, CD74
    '12': 'Memory/Effector T cells', # GZMK, CD27, CD3D, CD3E, TCF7, CCR7
    '13': 'Erythrocytes'            # HBA1, HBB, HBA2, SLC25A37, ALAS2
}

adata_proc.obs['cell_type'] = adata_proc.obs['leiden_0.3'].map(cluster_map)

# Check value counts
print(adata_proc.obs['cell_type'].value_counts())

# Generate labels.csv
labels_df = pd.DataFrame({
    'barcode': adata_proc.obs_names,
    'cell_type': adata_proc.obs['cell_type']
})

labels_df.to_csv('labels.csv', index=False)
print("Saved labels.csv with shape:", labels_df.shape)
print("Head:\n", labels_df.head(10))
print("Missing values in labels.csv:", labels_df.isnull().sum().to_dict())