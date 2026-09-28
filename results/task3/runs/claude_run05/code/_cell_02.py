import scanpy as sc
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
print(adata)
print(adata.X[:5,:5])
print(adata.obs.head())
print(adata.var.head())