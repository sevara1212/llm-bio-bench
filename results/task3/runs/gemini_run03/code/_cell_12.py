# Check CD14 vs FCGR3A in cluster 2 vs 4
for gene in ['CD14', 'FCGR3A', 'MS4A7', 'S100A8', 'S100A9']:
    print(gene, "in raw var?", gene in adata.raw.var_names)
    if gene in adata.raw.var_names:
        c2_mean = adata.raw[adata.obs['leiden'] == '2', gene].X.mean()
        c4_mean = adata.raw[adata.obs['leiden'] == '4', gene].X.mean()
        print(f"{gene}: Cl 2 = {c2_mean:.3f}, Cl 4 = {c4_mean:.3f}")