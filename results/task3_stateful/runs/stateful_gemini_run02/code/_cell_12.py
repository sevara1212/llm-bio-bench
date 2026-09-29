# Let's inspect CD4 vs CD8 markers across clusters at res=1.0:
for gene in ['CD3D', 'CD4', 'IL7R', 'CD8A', 'CD8B', 'NKG7', 'GNLY', 'MS4A1', 'CD14', 'FCGR3A', 'FCER1A', 'PPBP']:
    vals = []
    for c in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
        mean_val = adata[adata.obs['leiden'] == c, gene].X.mean()
        vals.append(f"{c}:{mean_val:.2f}")
    print(f"{gene:8s}: " + " ".join(vals))