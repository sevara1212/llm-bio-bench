# Let's see all 19 clusters' top 10 markers and mean expression of key canonical PBMC markers!
key_genes = ['CD3D', 'CD4', 'IL7R', 'CD8A', 'GNLY', 'NKG7', 'MS4A1', 'CD79A', 'CD14', 'FCGR3A', 'FCER1A', 'PPBP', 'HBB', 'LILRA4', 'FOXP3', 'GZMK', 'CCR7']
means = []
for c in sorted(adata_proc.obs['leiden_0.5'].unique(), key=lambda x: int(x)):
    cells = adata_proc.obs['leiden_0.5'] == c
    sub = adata_proc[cells, key_genes].X
    if hasattr(sub, 'toarray'):
        sub = sub.toarray()
    means.append(sub.mean(axis=0))

df_means = pd.DataFrame(means, index=[f"cl_{i}" for i in range(len(means))], columns=key_genes)
print(df_means.round(2))

# Also print top 10 genes for all clusters
for i in range(19):
    genes = adata_proc.uns['rank_genes_groups']['names'][str(i)][:7]
    print(f"Cluster {i}: {', '.join(genes)}")