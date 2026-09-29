import scanpy as sc

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

# Let's check clusters 0, 4, 6, 8, 11, 12, 14, 15, 16, 17
clusters_to_inspect = [0, 4, 6, 8, 11, 12, 14, 15, 16, 17]
genes = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'IL7R', 'GZMK', 'GZMB', 'GZMA', 'NKG7', 'NCAM1', 'FCGR3A', 'GNLY', 'TRAC', 'TRDC', 'FOXP3', 'IL2RA', 'KLRC1', 'XCL1', 'CLEC9A', 'CD1C', 'CLEC10A', 'CD34', 'SOX4', 'MKI67', 'TOP2A', 'GATA3', 'IL1RL1', 'FCER1A', 'HDC']

import pandas as pd
res_df = pd.DataFrame(index=[f"C{c}" for c in clusters_to_inspect])
for g in genes:
    if g in adata.var_names:
        vals = [float(adata[adata.obs['leiden'] == str(c), g].X.mean()) for c in clusters_to_inspect]
        res_df[g] = vals

pd.set_option('display.max_columns', 30)
pd.set_option('display.width', 1000)
print(res_df.round(2))