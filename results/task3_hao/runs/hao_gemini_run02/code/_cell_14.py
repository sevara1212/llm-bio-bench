# Let's inspect CD4, CD8A, CD8B, CCR7, S100A4, FOXP3, FCGR3A, CD14, HLA-DRA, CLEC10A, FCER1A across all clusters
for gene in ['CD4', 'CD8A', 'CD8B', 'CCR7', 'S100A4', 'FOXP3', 'FCGR3A', 'CD14', 'FCER1A', 'CLEC10A', 'CLEC9A', 'LILRA4', 'NKG7', 'GNLY', 'NCAM1', 'TRDC', 'TRGC1', 'TRGC2']:
    if gene in adata.var_names:
        vals = [f"{cl}:{np.mean(adata[adata.obs['leiden'] == cl, gene].X):.2f}" for cl in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x))]
        print(f"{gene:8s}: {' '.join(vals)}")