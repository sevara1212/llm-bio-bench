# Let's check non-ribosomal markers for cluster 0 and cluster 11
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')
for cl in ['0', '2', '3', '6', '16', '17', '11']:
    names = adata.uns['rank_genes_groups']['names'][cl]
    non_rp = [g for g in names if not g.startswith(('RP', 'MT-'))][:10]
    print(f"Cluster {cl} non-RP markers: {non_rp}")

# Also check CD4 vs CD8 vs Treg / MAIT / gamma-delta
check_genes = ['CD3D', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'IL7R', 'S100A4', 'FOXP3', 'KLRB1', 'SLC4A10', 'TRDV2', 'TRGV9', 'GZMK', 'GZMB', 'GNLY', 'NKG7', 'PRF1']
for cl in ['0', '2', '3', '6', '16', '17']:
    vals = {g: round(float(adata.raw[adata.obs['leiden_0.5'] == cl, g].X.mean()), 3) for g in check_genes if g in adata.raw.var_names}
    print(f"Cluster {cl}: {vals}")