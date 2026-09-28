import scanpy as sc
adata = sc.read_h5ad('clustered.h5ad')

# use raw (log-normalized) data for marker gene detection
adata_raw = adata.raw.to_adata()
adata_raw.obs['leiden'] = adata.obs['leiden']

sc.tl.rank_genes_groups(adata_raw, 'leiden', method='wilcoxon')

for cl in adata_raw.obs['leiden'].cat.categories:
    names = adata_raw.uns['rank_genes_groups']['names'][cl][:15]
    print(cl, list(names))