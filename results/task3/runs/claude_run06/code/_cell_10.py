import scanpy as sc
adata = sc.read_h5ad('clustered.h5ad')

# compare cluster 0 vs 6 directly
sub = adata[adata.obs['leiden'].isin(['0','6'])].copy()
sc.tl.rank_genes_groups(sub, 'leiden', groups=['0'], reference='6', method='wilcoxon')
res = sub.uns['rank_genes_groups']
print("Cluster0 vs 6 top genes (0 up):")
for i in range(20):
    print(res['names']['0'][i], round(res['logfoldchanges']['0'][i],2))