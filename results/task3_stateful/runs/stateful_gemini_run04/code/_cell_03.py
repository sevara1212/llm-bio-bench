adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=True)
print(adata)
print("Obs sample:", adata.obs_names[:5])
print("Var sample:", adata.var_names[:5])