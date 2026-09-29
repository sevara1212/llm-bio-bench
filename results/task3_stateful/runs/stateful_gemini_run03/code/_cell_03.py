adata = sc.read_10x_mtx('./filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
print(adata)
print("Barcodes sample:", adata.obs_names[:5].tolist())
print("Genes sample:", adata.var_names[:5].tolist())