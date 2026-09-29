# Let's check marker genes for standard PBMC types
# CD4 T cells: CD3D, CD3E, CD4, IL7R
# CD8 T cells: CD3D, CD8A, CD8B
# NK cells: GNLY, NKG7, NCAM1 (CD56)
# B cells: MS4A1 (CD20), CD19, CD79A
# CD14+ Monocytes: CD14, LYZ, S100A9, VCAN
# FCGR3A+ / CD16+ Monocytes: FCGR3A, MS4A7
# Dendritic cells (cDC): FCER1A, CST3, HLA-DQA1
# pDC: LILRA4, CLEC4C, TCF4, IL3RA
# Platelets / Megakaryocytes: PPBP, PF4
# Plasma cells: MZB1, SDC1, JCHAIN

# Let's compute rank_genes_groups on raw or log1p normalized data
adata_proc.obs['leiden'] = adata_hvg.obs['leiden']
sc.tl.rank_genes_groups(adata_proc, groupby='leiden', method='wilcoxon')

result = adata_proc.uns['rank_genes_groups']
groups = result['names'].dtype.names

for group in groups:
    top_genes = [result['names'][group][i] for i in range(10)]
    print(f"Cluster {group}: {', '.join(top_genes)}")