import scanpy as sc
adata = sc.read_h5ad('clustered.h5ad')

sc.tl.rank_genes_groups(adata, groupby='leiden_05', method='wilcoxon', use_raw=True)

for cl in adata.obs['leiden_05'].cat.categories:
    genes = adata.uns['rank_genes_groups']['names'][cl][:15]
    print(cl, list(genes))