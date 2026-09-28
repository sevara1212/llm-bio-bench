import scanpy as sc
import numpy as np
import pandas as pd

raw_barcodes = pd.read_csv('filtered_gene_bc_matrices/hg19/barcodes.tsv', header=None)[0].values

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Standard filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[(adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5), :].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=0)

sc.tl.leiden(adata, resolution=0.8, random_state=0)

genes_to_check = ['CD14', 'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'CD3D', 'CD4', 'CD8A', 'GNLY', 'NKG7', 'MS4A1', 'PPBP']
expr = pd.DataFrame(adata.raw[:, genes_to_check].X.toarray(), columns=genes_to_check, index=adata.obs_names)
expr['cluster'] = adata.obs['leiden']
print("\nMean expression by cluster (res 0.8):")
print(expr.groupby('cluster').mean().round(3))

# Check cluster names
# 0: T cells (CD3D high)
# 1: NK / Cytotoxic T cells (NKG7 high, GNLY high)
# 2: CD14+ Monocytes (CD14 high, LYZ, S100A9)
# 3: B cells (MS4A1 high, CD79A high)
# 4: FCGR3A+ / Non-classical Monocytes (FCGR3A / MS4A7 / LST1 high)
# 5: Dendritic cells (FCER1A, HLA-DRA high)
# 6: Megakaryocytes / Platelets (PPBP, PF4 high)