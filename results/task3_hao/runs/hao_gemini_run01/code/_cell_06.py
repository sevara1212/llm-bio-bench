import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("Total cells:", adata.n_obs)
print("Genes < 200:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("Genes < 300:", (adata.obs['n_genes_by_counts'] < 300).sum())
print("Genes < 500:", (adata.obs['n_genes_by_counts'] < 500).sum())
print("Genes > 5000:", (adata.obs['n_genes_by_counts'] > 5000).sum())
print("Counts > 40000:", (adata.obs['total_counts'] > 40000).sum())
print("pct_counts_mt > 10%:", (adata.obs['pct_counts_mt'] > 10).sum())
print("pct_counts_mt > 15%:", (adata.obs['pct_counts_mt'] > 15).sum())

# Gene detection frequency
n_cells_per_gene = np.array((adata.X > 0).sum(axis=0)).flatten()
print("Genes expressed in < 3 cells:", np.sum(n_cells_per_gene < 3))
print("Genes expressed in >= 3 cells:", np.sum(n_cells_per_gene >= 3))