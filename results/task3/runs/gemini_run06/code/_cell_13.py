# Let's inspect cluster 0 vs 6 and cluster 3 vs 5
# Check CD4, CD8A, CD8B, IL7R, GNLY, NKG7, FCGR3A, MS4A7, etc.
raw_adata = adata.raw.to_adata()
for gene in ['CD3D', 'CD3E', 'IL7R', 'CD4', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'CD14', 'FCGR3A', 'MS4A1', 'PPBP', 'FCER1A']:
    if gene in raw_adata.var_names:
        print(f"\nGene: {gene}")
        for cl in sorted(adata.obs['leiden_1'].unique(), key=lambda x: int(x)):
            idx = adata.obs['leiden_1'] == cl
            val = np.expm1(raw_adata[idx, gene].X.toarray()).mean() # or log mean
            logval = raw_adata[idx, gene].X.toarray().mean()
            print(f"  Clust {cl}: {logval:.3f}", end="")
        print()