import scanpy as sc
import numpy as np
import pandas as pd

# Load raw counts
adata = sc.read_h5ad('raw_counts.h5ad')
adata.layers['counts'] = adata.X.copy()

# QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Normalization & log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# PCA
sc.tl.pca(adata, mask_var="highly_variable", svd_solver='arpack', random_state=42)

# Neighborhood graph and UMAP
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata, random_state=42)

# Leiden clustering
sc.tl.leiden(adata, resolution=0.5, random_state=42, key_added='leiden')

# Check expression of key lineage markers
lineage_genes = ['CD3D', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'S100A4', 'FOXP3', 'NKG7', 'GNLY', 
                 'MS4A1', 'CD14', 'FCGR3A', 'CST3', 'FCER1A', 'CLEC9A', 'LILRA4', 'PPBP', 'MZB1', 'HBB', 'CD34', 'MKI67']
for g in lineage_genes:
    vals = [f"{cl}:{np.mean(adata[adata.obs['leiden'] == str(cl), g].X):.2f}" for cl in range(19)]
    print(f"{g:8s}: {' '.join(vals)}")

# Cluster annotations based on marker analysis:
# Cluster 0: CD4+ Naive T cells (CD3D+, CD4+, CCR7+, ribosomal high)
# Cluster 1: B cells (MS4A1, CD79A, CD79B, BANK1)
# Cluster 2: CD4+ Memory T cells (CD3D+, CD4+, IL7R+, S100A4+)
# Cluster 3: CD8+ T cells (CD3D+, CD8A+, CD8B+, CCL5, NKG7, CST7)
# Cluster 4: Proliferating cells (MKI67+, STMN1+, TOP2A+, TYMS+)
# Cluster 5: CD4+ T cells / Memory T cells (IL7R+, KLRB1+, CD69+)
# Cluster 6: Classical Monocytes (CD14+, CST3+, HLA-DR+, LYZ+)
# Cluster 7: Non-classical Monocytes (FCGR3A+, MS4A7+, LST1+, AIF1+)
# Cluster 8: NK cells (GNLY+, NKG7+, KLRD1+, CD3D-)
# Cluster 9: Plasmacytoid Dendritic Cells (pDC) (LILRA4+, TCF4+, IL3RA+, CCDC50+)
# Cluster 10: Plasma cells (MZB1+, JCHAIN+, SDC1+)
# Cluster 11: Hematopoietic stem/progenitor cells (HSPCs) (CD34+, CDK6+, SOX4+)
# Cluster 12: CD14+ Monocytes (CD14+, S100A9+, S100A8+, VCAN+, FCN1+)
# Cluster 13: Platelets (PPBP+, PF4+, TUBB1+)
# Cluster 14: NK cells / Cytotoxic cells (FCGR3A+, PRF1+, NKG7+, GNLY+)
# Cluster 15: Conventional Dendritic Cells (cDC) (CLEC9A+, CD74+, HLA-DRA+, CST3+)
# Cluster 16: T cells / MAIT / Treg / gdT cells (CD3+, KLRB1+, TNFRSF18+, IL7R+)
# Cluster 17: CD8+ Effector Memory T cells (CD3D+, CD8+, GZMK+, CD27+)
# Cluster 18: Erythrocytes (HBB+, HBA1+, HBA2+, ALAS2+)

cluster_labels = {
    '0': 'CD4+ Naive T cells',
    '1': 'B cells',
    '2': 'CD4+ Memory T cells',
    '3': 'CD8+ T cells',
    '4': 'Proliferating cells',
    '5': 'CD4+ T cells',
    '6': 'CD14+ Monocytes',
    '7': 'FCGR3A+ Monocytes',
    '8': 'NK cells',
    '9': 'Plasmacytoid Dendritic Cells',
    '10': 'Plasma cells',
    '11': 'HSPCs',
    '12': 'CD14+ Monocytes',
    '13': 'Platelets',
    '14': 'NK cells',
    '15': 'Conventional Dendritic Cells',
    '16': 'T cells',
    '17': 'CD8+ T cells',
    '18': 'Erythrocytes'
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_labels)

df_out = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type']
})

df_out.to_csv('labels.csv', index=False)
print(f"Successfully saved labels.csv with shape {df_out.shape}")
print(df_out['cell_type'].value_counts())