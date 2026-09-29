# Let's check marker gene expression across clusters
# Key PBMC markers:
# CD4+ T cells: CD3D, CD3E, CD4, IL7R
# CD8+ T cells: CD8A, CD8B
# NK cells: GNLY, NKG7, NCAM1 (CD56)
# B cells: MS4A1 (CD20), CD19, CD79A
# Monocytes (CD14+): CD14, LYZ
# Monocytes (FCGR3A+ / CD16+): FCGR3A, MS4A7
# Dendritic cells: FCER1A, CST3, HLA-DQA1, CLEC10A
# Plasmacytoid DC: IL3RA, CLEC4C, LILRA4
# Platelets / Megakaryocytes: PPBP, PF4
# Plasma cells: SDC1, MZB1, JCHAIN

sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

for g in groups:
    top_genes = [result['names'][g][i] for i in range(10)]
    print(f"Cluster {g}: {', '.join(top_genes)}")