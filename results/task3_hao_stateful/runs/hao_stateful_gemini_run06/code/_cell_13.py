# Let's inspect canonical markers across clusters
canonical = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'NCAM1', 'FCGR3A', 'NCAM1', 'NKG7', 'GNLY', 
             'CD14', 'FCGR3A', 'MS4A1', 'CD19', 'CD79A', 'PPBP', 'PF4', 'MZB1', 'JCHAIN', 'TCF4', 
             'IL3RA', 'CLEC9A', 'CD1C', 'HBB', 'MKI67']

# Let's check expression of CD4, CD8A, CD8B, CD14, FCGR3A, FCER1A, CLEC10A, CLEC9A across clusters
sc.tl.dendrogram(adata_proc, groupby='leiden_0.5')
check_genes = ['CD3D', 'CD4', 'CD8A', 'CD8B', 'NCAM1', 'FCGR3A', 'NKG7', 'GNLY', 'CD14', 'MS4A1', 'MZB1', 'TCF4', 'IL3RA', 'CLEC9A', 'CD1C', 'PPBP', 'HBB', 'MKI67']
for g in check_genes:
    if g in adata_proc.var_names:
        expr = [adata_proc[adata_proc.obs['leiden_0.5'] == str(i), g].X.mean() for i in range(19)]
        print(f"{g:8s}: " + " ".join([f"{v:4.2f}" for v in expr]))