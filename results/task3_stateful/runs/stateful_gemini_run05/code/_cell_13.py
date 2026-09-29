# Let's inspect CD8A, GNLY, NKG7, FCGR3A, MS4A1, LYZ, IL7R for each cluster
for g in groups:
    print(f"\n--- Cluster {g} ---")
    for gene in ['CD3D', 'CD4', 'IL7R', 'CD8A', 'GNLY', 'NKG7', 'MS4A1', 'CD14', 'FCGR3A', 'FCER1A', 'PPBP']:
        if gene in adata_filtered.raw.var_names:
            idx = adata_filtered.raw.var_names.get_loc(gene)
            val = np.mean(adata_filtered.raw.X[(adata_filtered.obs['leiden'] == g).values, idx].toarray())
            print(f"  {gene}: {val:.2f}", end="")
    print()