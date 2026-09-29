# Let's check marker genes for standard PBMC types
marker_genes = ['IL7R', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'S100A4', 'MS4A1', 'CD79A', 'CD14', 'LYZ', 'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP']

# Check which markers are in adata
print([m for m in marker_genes if m in adata.var_names])

# Compute rank genes groups
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)