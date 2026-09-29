marker_genes = ['IL7R', 'CD4', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 'GNLY', 'NKG7', 'FCGR3A', 'MS4A7', 'CST3', 'PPBP', 'FCER1A', 'HLA-DQA1']
# Let's inspect top 10 genes for each cluster
for g in groups:
    top_genes = [result['names'][g][i] for i in range(10)]
    print(f"Cluster {g}: {top_genes}")