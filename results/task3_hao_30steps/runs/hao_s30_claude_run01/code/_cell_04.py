import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, percent_top=None)

for q in [0.9,0.95,0.98,0.99,0.995,0.999]:
    print(q, 'n_genes', np.quantile(adata.obs.n_genes_by_counts,q), 'total_counts', np.quantile(adata.obs.total_counts,q))

print('min genes', adata.obs.n_genes_by_counts.min())
print('n cells with pct_mt>10', (adata.obs.pct_counts_mt>10).sum())
print('n cells with pct_mt>15', (adata.obs.pct_counts_mt>15).sum())