# Let's inspect gene filtering and cell filtering
import scanpy as sc
adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
sc.pp.filter_genes(adata, min_cells=3)
print("After filter_genes min_cells=3:", adata)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Look at distribution of pct_counts_mt
print("mt > 15%:", (adata.obs['pct_counts_mt'] > 15).sum())
print("mt > 10%:", (adata.obs['pct_counts_mt'] > 10).sum())
print("n_genes < 200:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("n_genes < 500:", (adata.obs['n_genes_by_counts'] < 500).sum())