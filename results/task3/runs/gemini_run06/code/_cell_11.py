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

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_sub = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_sub, max_value=10)
sc.pp.pca(adata_sub, svd_solver='arpack', random_state=0)
sc.pp.neighbors(adata_sub, n_neighbors=10, n_pcs=40, random_state=0)

# Check resolution 0.8
sc.tl.leiden(adata_sub, resolution=0.8, key_added='leiden_08', random_state=0)
adata.obs['leiden_08'] = adata_sub.obs['leiden_08']
sc.tl.rank_genes_groups(adata, 'leiden_08', method='wilcoxon')
result = adata.uns['rank_genes_groups']
print("Markers res 0.8:")
for g in result['names'].dtype.names:
    print(f"Cluster {g} (n={sum(adata.obs['leiden_08']==g)}): {list(result['names'][g][:6])}")

# Check key marker gene expression means
marker_genes = ['IL7R', 'CD3D', 'CD14', 'LYZ', 'MS4A1', 'CD79A', 'GNLY', 'NKG7', 'FCGR3A', 'MS4A7', 'CST3', 'FCER1A', 'PPBP', 'CD8A']
print("\nMean expression res 0.8:")
mean_df = pd.DataFrame(index=marker_genes)
for g in sorted(adata.obs['leiden_08'].unique()):
    cells = adata[adata.obs['leiden_08'] == g]
    mean_df[f'Clust_{g}'] = [cells[:, gene].raw.X.mean() if gene in adata.raw.var_names else 0 for gene in marker_genes]
print(mean_df.round(2))