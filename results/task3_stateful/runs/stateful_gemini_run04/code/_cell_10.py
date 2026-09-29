# Check genes like CD4, CD8A, CD8B, IL7R, CCR7, S100A4, GNLY, NKG7, MS4A1, CD14, FCGR3A, FCER1A, PPBP
test_genes = ['IL7R', 'CCR7', 'S100A4', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP']
for g in test_genes:
    print(g, g in adata_filtered.raw.var_names)