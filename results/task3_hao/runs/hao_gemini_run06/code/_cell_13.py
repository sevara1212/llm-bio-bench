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
sc.tl.leiden(adata, resolution=0.5, random_state=42, flavor='igraph')
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

res = adata.uns['rank_genes_groups']
for c in range(19):
    c_str = str(c)
    top_genes = [res['names'][c_str][i] for i in range(7)]
    sub = adata[adata.obs['leiden'] == c_str]
    c_cd4 = sub[:, 'CD4'].X.mean()
    c_cd8a = sub[:, 'CD8A'].X.mean()
    c_cd8b = sub[:, 'CD8B'].X.mean() if 'CD8B' in adata.var_names else 0
    c_cd3d = sub[:, 'CD3D'].X.mean()
    c_cd14 = sub[:, 'CD14'].X.mean()
    c_fcgr3a = sub[:, 'FCGR3A'].X.mean()
    c_nkg7 = sub[:, 'NKG7'].X.mean()
    c_ms4a1 = sub[:, 'MS4A1'].X.mean()
    print(f"C{c:2d} (n={sub.shape[0]}): top={top_genes} | CD3D={c_cd3d:.1f}, CD4={c_cd4:.1f}, CD8A={c_cd8a:.1f}, CD8B={c_cd8b:.1f}, MS4A1={c_ms4a1:.1f}, CD14={c_cd14:.1f}, FCGR3A={c_fcgr3a:.1f}, NKG7={c_nkg7:.1f}")