# Let's inspect standard Seurat / Scanpy PBMC3k tutorial filtering and clustering
# In Scanpy pbmc3k tutorial:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# adata.raw = adata
# adata = adata[:, adata.var.highly_variable]
# sc.pp.scale(adata, max_value=10)
# sc.pp.pca(adata, svd_solver='arpack')
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40) # or standard
# sc.tl.leiden(adata)
import scanpy as sc

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print("After gene filter min_cells=3:", adata.shape)

adata_filtered = adata[(adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5), :].copy()
print("After QC filter:", adata_filtered.shape)
print("Barcodes retained:", adata_filtered.n_obs)