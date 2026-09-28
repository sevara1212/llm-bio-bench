import scanpy as sc
import matplotlib
matplotlib.use('Agg')

adata = sc.read_h5ad('raw.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("Before filtering:", adata.shape)

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print("After min_genes/min_cells filter:", adata.shape)

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

print("After QC filtering:", adata.shape)

adata.write('qc.h5ad')