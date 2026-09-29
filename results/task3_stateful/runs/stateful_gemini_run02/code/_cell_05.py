# Let's see how many cells are filtered with standard cutoffs
# In Seurat pbmc3k tutorial:
# subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
print("n_genes < 200:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("n_genes > 2500:", (adata.obs['n_genes_by_counts'] > 2500).sum())
print("pct_counts_mt > 5:", (adata.obs['pct_counts_mt'] > 5).sum())

# In scanpy pbmc3k tutorial:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]