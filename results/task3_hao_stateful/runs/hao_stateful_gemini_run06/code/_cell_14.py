# Let's inspect cluster 4, 11, 14, 16 more closely:
genes = ['CD4', 'CD8A', 'CD8B', 'TRAC', 'TRDC', 'TRGC1', 'FOXP3', 'IL2RA', 'CCR7', 'SELL', 'IL7R', 'GZMK', 'GZMB', 'TYMS', 'MKI67', 'TCL1A']
for g in genes:
    if g in adata_proc.var_names:
        expr = [adata_proc[adata_proc.obs['leiden_0.5'] == str(i), g].X.mean() for i in [0, 1, 2, 4, 11, 14, 16, 18]]
        print(f"{g:8s}: " + " ".join([f"{v:4.2f}" for v in expr]))

print("\nCluster 11 top markers:")
print(list(markers['11']))
print("\nCluster 14 top markers:")
print(list(markers['14']))
print("\nCluster 4 top markers:")
print(list(markers['4']))