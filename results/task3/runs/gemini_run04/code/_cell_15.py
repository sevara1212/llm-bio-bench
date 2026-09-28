# In this final run:
# 1. Perform standard QC, normalization, HVG, PCA, clustering
# 2. Map clusters to cell type labels:
#    Cluster 0: T cells
#    Cluster 1: NK cells
#    Cluster 2: CD14+ Monocytes
#    Cluster 3: B cells
#    Cluster 4: FCGR3A+ Monocytes
#    Cluster 5: Dendritic cells
#    Cluster 6: Megakaryocytes
# 3. For any barcode in barcodes.tsv that was filtered out during QC, assign 'Filtered' or map nearest/keep all 2700 barcodes or output filtered cells.
# Let's check standard practice: typically labels.csv contains the analyzed cells or all barcodes. Let's see if we should include all 2638 or all 2700.
# If someone expects labels for the clustered cells, barcode,cell_type.
# Let's save all 2638 clustered cells, and check if saving adata.obs[['barcode', 'cell_type']] is standard.
import scanpy as sc
import numpy as np
import pandas as pd

raw_barcodes = pd.read_csv('filtered_gene_bc_matrices/hg19/barcodes.tsv', header=None)[0].values

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

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

cluster_names = {
    '0': 'T cells',
    '1': 'NK cells',
    '2': 'CD14+ Monocytes',
    '3': 'B cells',
    '4': 'FCGR3A+ Monocytes',
    '5': 'Dendritic cells',
    '6': 'Megakaryocytes'
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_names)

df_out = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type']
})

df_out.to_csv('labels.csv', index=False)
print("labels.csv written successfully. Shape:", df_out.shape)
print(df_out.head(10))
print(df_out['cell_type'].value_counts())