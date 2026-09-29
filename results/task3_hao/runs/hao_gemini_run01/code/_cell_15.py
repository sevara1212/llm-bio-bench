# Let's inspect cluster 12 and 17, and create the final cell type mapping and save labels.csv
import scanpy as sc
import pandas as pd
import numpy as np

# We have 1 execution left, let's write a complete, self-contained, robust script:
# 1. Load data
# 2. QC & preprocessing
# 3. Clustering
# 4. Marker identification & annotation
# 5. Output labels.csv
# 6. Verify labels.csv

adata = sc.read_h5ad('raw_counts.h5ad')
barcodes = adata.obs_names.copy()

# Filter genes
sc.pp.filter_genes(adata, min_cells=3)

# Normalize and log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# HVG, PCA, Neighbors, Leiden
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata_hvg, resolution=0.5, random_state=42)
adata.obs['leiden'] = adata_hvg.obs['leiden']

# Let's check marker genes for each cluster
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
result = adata.uns['rank_genes_groups']

# Cluster annotations based on markers:
# Cluster 0: CD4+ Naive T cells (LEF1, TCF7, CD3D, CCR7, MAL, IL7R) -> CD4+ T cell
# Cluster 1: B cells (MS4A1, BANK1, CD79A, RALGPS2, CD79B, CD37) -> B cell
# Cluster 2: CD8+ Effector/Cytotoxic T cells (CCL5, NKG7, GZMH, CST7, CD3D, CD8A, PRF1) -> CD8+ T cell
# Cluster 3: CD4+ Memory T cells (IL7R, IL32, LTB, TRAC, CD3D, CD69, CD4) -> CD4+ T cell
# Cluster 4: Proliferating cells (STMN1, HMGB2, DEK, TYMS, PCNA, MKI67) -> Proliferating T cell
# Cluster 5: CD14+ Monocytes (CD14, VCAN, FTL, S100A8, S100A9, MNDA) -> CD14+ Monocyte
# Cluster 6: MAIT / CD8+ Memory T cells (KLRB1, GZMK, DUSP2, SLC4A10, NCR3, CD3D) -> MAIT cell
# Cluster 7: FCGR3A+ / Non-classical Monocytes (FCGR3A, MS4A7, LST1, AIF1, CFD) -> FCGR3A+ Monocyte
# Cluster 8: NK cells (GNLY, KLRD1, XCL1, XCL2, KLRF1, KLRC1, CTSW, no CD3D) -> NK cell
# Cluster 9: Conventional Dendritic Cells / cDC2 (HLA-DRA, HLA-DQB1, CD1C, FCER1A, LYZ, CST3) -> Dendritic cell
# Cluster 10: Plasmacytoid Dendritic Cells / pDC (TCF4, CCDC50, PLD4, IL3RA, LILRA4, IRF7) -> Plasmacytoid dendritic cell
# Cluster 11: Plasma cells (MZB1, JCHAIN, TNFRSF17, SDC1, HSP90B1) -> Plasma cell
# Cluster 12: Progenitor / HSPC / Pre-B cells (SOX4, CDK6, ZFAS1, PRSS57, ANKRD28) -> Hematopoietic stem and progenitor cell
# Cluster 13: Platelets / Megakaryocytes (PF4, PPBP, CAVIN2, NRGN, GNG11, TUBB1) -> Platelet
# Cluster 14: Cytotoxic NK cells (PRF1, NKG7, KLRF1, CST7, KLRD1, FGFBP2, GZMB, GNLY) -> NK cell
# Cluster 15: CD4+ T cells / Memory CD4+ T (CD3D, GZMK, CD27, TCF7, LCK) -> CD4+ T cell
# Cluster 16: Conventional Dendritic Cells / cDC1 (CLEC9A, WDFY4, CPNE3, SNX3, HLA-DPA1) -> Conventional dendritic cell
# Cluster 17: gamma-delta T / Activated T cells (TRDC, TNFRSF18, TNFRSF4, KLRB1, CD7) -> gamma-delta T cell
# Cluster 18: Plasmacytoid dendritic cells / pDC subset (LILRA4, PLD4, PPP1R14A, SAMHD1, APP) -> Plasmacytoid dendritic cell

cluster_annotation = {
    '0': 'CD4+ Naive T cell',
    '1': 'B cell',
    '2': 'CD8+ Cytotoxic T cell',
    '3': 'CD4+ Memory T cell',
    '4': 'Proliferating T cell',
    '5': 'CD14+ Monocyte',
    '6': 'MAIT cell',
    '7': 'FCGR3A+ Monocyte',
    '8': 'NK cell',
    '9': 'Conventional dendritic cell',
    '10': 'Plasmacytoid dendritic cell',
    '11': 'Plasma cell',
    '12': 'Hematopoietic progenitor cell',
    '13': 'Platelet',
    '14': 'NK cell',
    '15': 'CD4+ T cell',
    '16': 'Conventional dendritic cell',
    '17': 'gamma-delta T cell',
    '18': 'Plasmacytoid dendritic cell'
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotation)

# Save labels.csv
df_out = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type']
})

df_out.to_csv('labels.csv', index=False)
print("labels.csv written successfully!")
print("Shape:", df_out.shape)
print("Distribution:\n", df_out['cell_type'].value_counts())
print("\nFirst 10 rows:\n", df_out.head(10))