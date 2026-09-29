# Let's inspect known PBMC markers across clusters:
# T cells: CD3D, CD3E, CD3G
# CD4 T cells: CD4, IL7R
# CD8 T cells: CD8A, CD8B
# NK cells: GNLY, NKG7, NCAM1 (CD56)
# B cells: MS4A1 (CD20), CD19, CD79A
# Monocytes: CD14, LYZ, FCGR3A (CD16), MS4A7
# Dendritic cells: FCER1A, CST3, CLEC4C, LILRA4
# Platelets / Megakaryocytes: PPBP, PF4
# Erythrocytes: HBB, HBA1, HBA2

markers = ['CD3D', 'CD3E', 'CD4', 'IL7R', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'MS4A1', 'CD19', 'CD79A', 'CD14', 'LYZ', 'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP', 'PF4', 'HBB']
present_markers = [m for m in markers if m in adata_proc.raw.var_names]

sc.tl.rank_genes_groups(adata_proc, groupby='leiden_0.5', method='wilcoxon')
top_genes = pd.DataFrame(adata_proc.uns['rank_genes_groups']['names']).head(5)
print("Top 5 marker genes per cluster (res 0.5):")
print(top_genes)