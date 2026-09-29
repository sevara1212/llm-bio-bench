# Let's inspect the markers per cluster without ribo/mito genes to see clear biological identity:
non_ribo_mt = ~adata_proc.var_names.str.startswith(('RPS', 'RPL', 'MT-'))
adata_sub = adata_proc[:, non_ribo_mt].copy()
sc.tl.rank_genes_groups(adata_sub, groupby='leiden_0.5', method='wilcoxon')
top_markers = pd.DataFrame(adata_sub.uns['rank_genes_groups']['names']).head(8)
print("Top markers (non-ribo/mt):")
for col in top_markers.columns:
    print(f"Cluster {col}: {top_markers[col].tolist()[:5]}")