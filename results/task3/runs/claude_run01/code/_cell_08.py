import scanpy as sc
adata = sc.read_h5ad('clustered.h5ad')

# use raw (log-normalized full gene set) for marker gene finding
adata_full = adata.raw.to_adata()
adata_full.obs['leiden'] = adata.obs['leiden']

sc.tl.rank_genes_groups(adata_full, 'leiden', method='wilcoxon')

for cl in adata_full.obs['leiden'].cat.categories:
    names = adata_full.uns['rank_genes_groups']['names'][cl][:15]
    print(cl, list(names))