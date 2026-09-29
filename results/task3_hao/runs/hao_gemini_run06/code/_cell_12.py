import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42)

for c in range(19):
    c_str = str(c)
    cells = adata.obs['leiden'] == c_str
    sub = adata[cells]
    
    # Check top markers
    sc.tl.rank_genes_groups(adata, 'leiden', groups=[c_str], reference='rest', method='wilcoxon')
    top_genes = [adata.uns['rank_genes_groups']['names'][c_str][i] for i in range(8)]
    
    # Check key lineage markers mean
    key_markers = ['CD3D', 'CD4', 'CD8A', 'IL7R', 'CCR7', 'GZMK', 'NKG7', 'GNLY', 'NCAM1', 'FCGR3A', 
                   'MS4A1', 'CD79A', 'MZB1', 'CD14', 'LYZ', 'S100A9', 'MS4A7', 'HLA-DRA', 'CLEC9A', 
                   'TCF4', 'IL3RA', 'PPBP', 'HBB', 'MKI67', 'CD34']
    mean_vals = {g: round(float(sub[:, g].X.mean()), 2) for g in key_markers if g in adata.var_names and sub[:, g].X.mean() > 0.5}
    print(f"Cluster {c} (n={cells.sum()}): Top: {', '.join(top_genes)} | High: {mean_vals}")