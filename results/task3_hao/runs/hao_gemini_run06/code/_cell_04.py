import scanpy as sc

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print("pct_counts_mt quantiles:")
print(adata.obs['pct_counts_mt'].quantile([0, 0.01, 0.05, 0.5, 0.95, 0.98, 0.99, 1.0]))
print("n_genes_by_counts quantiles:")
print(adata.obs['n_genes_by_counts'].quantile([0, 0.01, 0.05, 0.5, 0.95, 0.98, 0.99, 1.0]))
print("total_counts quantiles:")
print(adata.obs['total_counts'].quantile([0, 0.01, 0.05, 0.5, 0.95, 0.98, 0.99, 1.0]))

# Check if barcodes in labels.csv should include all 7841 barcodes or only filtered ones!
# "Save labels.csv with columns barcode,cell_type."