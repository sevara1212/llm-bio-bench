import scanpy as sc, pandas as pd, numpy as np
p='filtered_gene_bc_matrices/hg19'
adata=sc.read_10x_mtx(p,var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True,log1p=False)
print(adata)
print(adata.obs[['total_counts','n_genes_by_counts','pct_counts_mt']].describe(percentiles=[.01,.05,.1,.5,.9,.95,.99]).round(2))
print('mito high', (adata.obs.pct_counts_mt>10).sum(),(adata.obs.pct_counts_mt>15).sum(),(adata.obs.pct_counts_mt>20).sum())
# top gene
print(adata.var.sort_values('n_cells_by_counts',ascending=False)[['n_cells_by_counts','total_counts','mt']].head(15))