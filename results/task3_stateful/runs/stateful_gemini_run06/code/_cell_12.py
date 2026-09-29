# Let's run differential expression between cluster 0 and cluster 1
sc.tl.rank_genes_groups(adata_scaled, 'leiden', groups=['0'], reference='1', method='wilcoxon')
print("Top in 0 vs 1:", adata_scaled.uns['rank_genes_groups']['names']['0'][:10].tolist())
sc.tl.rank_genes_groups(adata_scaled, 'leiden', groups=['1'], reference='0', method='wilcoxon')
print("Top in 1 vs 0:", adata_scaled.uns['rank_genes_groups']['names']['1'][:10].tolist())

# Also check S100A4 and CCR7 distribution or IL7R
print("Mean IL7R, CCR7, S100A4, LDHB:")
for gene in ['IL7R', 'CCR7', 'S100A4', 'LDHB', 'MAL', 'TCF7', 'SELL', 'CD27', 'GZMA']:
    if gene in adata.raw.var_names:
        c0 = adata.raw[adata_scaled.obs['leiden'] == '0', gene].X.toarray().mean()
        c1 = adata.raw[adata_scaled.obs['leiden'] == '1', gene].X.toarray().mean()
        print(f"{gene:10s}: c0={c0:.3f}, c1={c1:.3f}")