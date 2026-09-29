# Let's inspect cluster markers for res 0.3 and check expression of canonical markers
sc.tl.rank_genes_groups(adata, groupby='leiden_0.3', method='wilcoxon')
markers = adata.uns['rank_genes_groups']
for group in sorted(markers['names'].dtype.names, key=int):
    top_genes = [markers['names'][group][i] for i in range(8)]
    print(f"Cluster {group}: {', '.join(top_genes)}")

canonical = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'IL7R', 'GNLY', 'NKG7', 'MS4A1', 'CD19', 'CD14', 'FCGR3A', 'CST3', 'FCER1A', 'PPBP', 'MZB1', 'LILRA4']
print("\nChecking canonical markers in adata.raw:")
for g in canonical:
    if g not in adata.raw.var_names:
        print(f"Missing: {g}")