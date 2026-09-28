sc.tl.rank_genes_groups(adata, 'leiden', groups=['0', '6'], reference='rest', method='wilcoxon')
print("Top in 0 vs rest:", adata.uns['rank_genes_groups']['names']['0'][:10])
print("Top in 6 vs rest:", adata.uns['rank_genes_groups']['names']['6'][:10])

# Check markers like CCR7, S100A4
for g in ['CCR7', 'S100A4', 'IL7R', 'MALAT1', 'FOXP3', 'CD27']:
    if g in adata.raw.var_names:
        e0 = adata.raw[:, g].X.toarray().flatten()[adata.obs['leiden'] == '0'].mean()
        e6 = adata.raw[:, g].X.toarray().flatten()[adata.obs['leiden'] == '6'].mean()
        print(f"{g}: cluster 0 = {e0:.3f}, cluster 6 = {e6:.3f}")