# Check the scanpy pbmc3k tutorial cluster names and mapping:
# In the Scanpy tutorial:
# new_cluster_names = [
#     'CD4 T', 'CD14+ Monocytes',
#     'B', 'CD8 T',
#     'FCGR3A+ Monocytes',
#     'NK', 'Dendritic',
#     'Megakaryocytes'
# ]
# adata.rename_categories('leiden', new_cluster_names)
# In scanpy pbmc3k tutorial, it had 8 clusters!
# Let's check what resolution gives exactly 8 clusters, or what resolution Scanpy tutorial uses:
# Actually in Scanpy tutorial: sc.tl.leiden(adata) default resolution is 1.0, but maybe with slightly different parameters it gave 8 clusters.
# Let's check resolution 0.8:
# Res 0.8 gave 7 clusters:
# 0: CD4 T cells / T cells (1181)
# 1: NK / CD8 T cells (429)
# 2: CD14+ Monocytes (417)
# 3: B cells (340)
# 4: FCGR3A+ Monocytes (222)
# 5: Dendritic cells (36)
# 6: Megakaryocytes (13)

# Let's write a script that inspects markers for res 1.0 and res 0.8
import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata_sub = adata.copy()
sc.pp.highly_variable_genes(adata_sub, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_sub = adata_sub[:, adata_sub.var.highly_variable].copy()
sc.pp.scale(adata_sub, max_value=10)
sc.pp.pca(adata_sub, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_sub, n_neighbors=10, n_pcs=40, random_state=0)
sc.tl.leiden(adata_sub, resolution=1.0, key_added='leiden_1', random_state=0)

adata.obs['leiden'] = adata_sub.obs['leiden_1']
raw_adata = adata.raw.to_adata()

genes = ['CD3D', 'CD3E', 'IL7R', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'CD14', 'LYZ', 'FCGR3A', 'MS4A7', 'MS4A1', 'CD79A', 'FCER1A', 'HLA-DRA', 'PPBP']
df_expr = pd.DataFrame(index=genes)
for cl in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
    idx = (adata.obs['leiden'] == cl).values
    df_expr[cl] = [raw_adata[idx, g].X.mean() for g in genes]
print(df_expr.round(3))