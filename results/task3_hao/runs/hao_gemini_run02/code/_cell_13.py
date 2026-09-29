# Let's check marker genes expression specifically for all clusters
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, mask_var="highly_variable", svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42, key_added='leiden')

check_genes = [
    'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'IL7R', 'S100A4', 'FOXP3',
    'NKG7', 'GNLY', 'NCAM1', 'FCGR3A', 'MS4A1', 'CD79A', 'CD14', 'LYZ', 'MS4A7',
    'FCER1A', 'CLEC9A', 'CLEC4C', 'IL3RA', 'LILRA4', 'TCF4', 'PPBP', 'PF4',
    'MZB1', 'JCHAIN', 'HBB', 'HBA1', 'MKI67', 'STMN1', 'CD34', 'CDK6'
]

# Calculate mean expression per cluster for check_genes
mean_df = pd.DataFrame(index=check_genes, columns=sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)))
for gene in check_genes:
    if gene in adata.var_names:
        for cl in mean_df.columns:
            mean_df.loc[gene, cl] = np.mean(adata[adata.obs['leiden'] == cl, gene].X)

print(mean_df.to_string())