# Let's inspect known PBMC markers across the clusters or find marker genes!
sc.tl.rank_genes_groups(adata_proc, groupby='leiden_0.5', method='wilcoxon')
markers = pd.DataFrame(adata_proc.uns['rank_genes_groups']['names']).head(10)
print(markers)

# Also check key PBMC canonical markers:
# CD4 T cells: CD3D, CD3E, CD4, IL7R
# CD8 T cells: CD8A, CD8B
# NK cells: GNLY, NKG7, NCAM1 (CD56), FCGR3A (CD16)
# B cells: MS4A1 (CD20), CD19, CD79A
# Monocytes: CD14 (Classical), FCGR3A / MS4A7 (Non-classical)
# Dendritic cells: FCER1A, CST3, HLA-DQA1, CLEC9A, LILRA4 / CLEC4C (pDC)
# Platelets / Megakaryocytes: PPBP, PF4
# Plasma cells: MZB1, SDC1 (CD138), JCHAIN
key_genes = ['CD3D', 'CD3E', 'CD4', 'IL7R', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'MS4A1', 'CD19', 'CD79A', 'CD14', 'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP', 'PF4', 'MZB1', 'LILRA4']
present_keys = [g for g in key_genes if g in adata_proc.var_names]
print("\nPresent key genes:", present_keys)