import scanpy as sc
adata = sc.read_h5ad('adata_ranked.h5ad')

# Compare cluster 0 vs 6 directly
sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='6', method='wilcoxon', key_added='0v6')
res = adata.uns['0v6']
names = res['names']['0'][:20]
lfc = res['logfoldchanges']['0'][:20]
print("Cluster0 up vs 6:", list(zip(names, lfc.round(2))))

sc.tl.rank_genes_groups(adata, 'leiden', groups=['6'], reference='0', method='wilcoxon', key_added='6v0')
res2 = adata.uns['6v0']
names2 = res2['names']['6'][:20]
lfc2 = res2['logfoldchanges']['6'][:20]
print("Cluster6 up vs 0:", list(zip(names2, lfc2.round(2))))