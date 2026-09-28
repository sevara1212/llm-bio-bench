import scanpy as sc

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=False)
print("Initial shape:", adata.shape)

# Scanpy PBMC3k tutorial filtering:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

n_genes = (adata.X > 0).sum(axis=1)
print("min genes:", n_genes.min())
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print("cells with mt < 5% and n_genes < 2500:", ((adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5)).sum())