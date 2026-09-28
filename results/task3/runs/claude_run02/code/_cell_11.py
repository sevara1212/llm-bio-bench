import scanpy as sc
adata = sc.read_h5ad('adata_ranked.h5ad')

print(adata.obs.groupby('leiden')[['n_genes_by_counts','total_counts','pct_counts_mt']].mean())