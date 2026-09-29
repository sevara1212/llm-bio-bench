import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42)

# Check mean expression of key markers across clusters
markers = [
    'CD3D', 'CD4', 'CD8A', 'IL7R', 'CCR7', 'FOXP3',
    'NKG7', 'GNLY', 'NCAM1', 'FCGR3A',
    'MS4A1', 'CD79A', 'MZB1',
    'CD14', 'LYZ', 'S100A9', 'MS4A7',
    'HLA-DRA', 'CD1C', 'CLEC9A', 'TCF4', 'IL3RA', 'LILRA4',
    'PPBP', 'PF4',
    'HBB', 'HBA1',
    'MKI67', 'TOP2A', 'STMN1',
    'CD34', 'SOX4'
]

# Create a summary table
df_list = []
for c in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
    cells = adata.obs['leiden'] == c
    sub = adata[cells, :]
    means = {gene: float(sub[:, gene].X.mean()) if gene in adata.var_names else 0.0 for gene in markers}
    means['cluster'] = c
    means['n_cells'] = cells.sum()
    df_list.append(means)

df_summary = pd.DataFrame(df_list).set_index('cluster')
pd.set_option('display.max_columns', 35)
pd.set_option('display.width', 1000)
print(df_summary)