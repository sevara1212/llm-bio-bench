# Let's inspect cluster 11, cluster 4, cluster 5, cluster 14, etc. more closely
# Cluster 11: CDK6, ZFAS1, SOX4, ANKRD28, PRSS57, SERPINB1 - CD34? HSPC / Progenitor?
# Let's check stem/progenitor markers for cluster 11:
hspc_genes = ['CD34', 'SPINK2', 'AVP', 'SOX4', 'GATA2', 'MPO', 'PRSS57', 'CDK6']
print("Cluster 11 expression:")
for g in hspc_genes:
    if g in adata.var_names:
        print(f"  {g}: {adata[adata.obs['leiden'] == '11', g].X.mean():.3f} (all: {adata[:, g].X.mean():.3f})")

# Cluster 5: HLA-DRA, CD1C, FCER1A, CLEC10A? (cDC2 / dendritic cells)
dc_genes = ['CD1C', 'FCER1A', 'CLEC10A', 'HLA-DRA', 'CST3', 'ITGAX']
print("\nCluster 5 expression:")
for g in dc_genes:
    if g in adata.var_names:
        print(f"  {g}: {adata[adata.obs['leiden'] == '5', g].X.mean():.3f} (all: {adata[:, g].X.mean():.3f})")

# Cluster 15: CLEC9A, WDFY4, XCR1, CADM1 (cDC1)
print("\nCluster 15 expression:")
for g in ['CLEC9A', 'WDFY4', 'XCR1', 'CADM1', 'HLA-DRA', 'CST3']:
    if g in adata.var_names:
        print(f"  {g}: {adata[adata.obs['leiden'] == '15', g].X.mean():.3f} (all: {adata[:, g].X.mean():.3f})")

# Cluster 9: pDC markers: LILRA4, CLEC4C, TCF4, IRF7, IL3RA, PLD4
print("\nCluster 9 expression:")
for g in ['LILRA4', 'CLEC4C', 'TCF4', 'IRF7', 'IL3RA', 'PLD4']:
    if g in adata.var_names:
        print(f"  {g}: {adata[adata.obs['leiden'] == '9', g].X.mean():.3f} (all: {adata[:, g].X.mean():.3f})")

# Cluster 10: Plasma cells: SDC1 (CD138), MZB1, JCHAIN, TNFRSF17 (BCMA), PRDM1 (BLIMP1)
print("\nCluster 10 expression:")
for g in ['MZB1', 'JCHAIN', 'TNFRSF17', 'SDC1', 'PRDM1', 'CD38']:
    if g in adata.var_names:
        print(f"  {g}: {adata[adata.obs['leiden'] == '10', g].X.mean():.3f} (all: {adata[:, g].X.mean():.3f})")

# Cluster 14: T cells? CD4? Treg? FOXP3?
print("\nCluster 14 expression:")
for g in ['CD3D', 'CD4', 'CD8A', 'IL7R', 'FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF4', 'TNFRSF18', 'IKZF2']:
    if g in adata.var_names:
        print(f"  {g}: {adata[adata.obs['leiden'] == '14', g].X.mean():.3f} (all: {adata[:, g].X.mean():.3f})")