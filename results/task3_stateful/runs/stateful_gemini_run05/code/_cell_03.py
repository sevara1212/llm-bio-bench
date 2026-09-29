adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
print(adata)