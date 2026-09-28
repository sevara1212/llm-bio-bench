import scanpy as sc
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
print(adata)
print(adata.X[:10,:10])
print(adata.var.head())
print(adata.obs.head())