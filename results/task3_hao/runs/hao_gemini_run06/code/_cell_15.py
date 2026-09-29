# Let's inspect cluster 3 and 11, and write the labels.csv file carefully
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42, flavor='igraph')

# Inspect markers for C3 vs C11
sc.tl.rank_genes_groups(adata, 'leiden', groups=['3', '11', '16'], method='wilcoxon')
print("C3 top:", [adata.uns['rank_genes_groups']['names']['3'][i] for i in range(10)])
print("C11 top:", [adata.uns['rank_genes_groups']['names']['11'][i] for i in range(10)])
print("C16 top:", [adata.uns['rank_genes_groups']['names']['16'][i] for i in range(10)])

# Cell type mapping:
# C0: Proliferating T cells (MKI67+, TOP2A+, STMN1+, CD3D+)
# C1: Non-classical monocytes (FCGR3A/CD16+, LST1+, MS4A7+, AIF1+)
# C2: Platelets / Megakaryocytes (PPBP+, PF4+, CAVIN2+)
# C3: Conventional dendritic cells / cDC2 (HLA-DRA+, CST3+, CD74+, HLA-DP/DQ+)
# C4: CD8+ effector memory T cells / Memory CD8+ T cells (CD3D+, CD8A+, GZMK+, IL7R+, KLRB1+)
# C5: Erythrocytes (HBB+, HBA1+, HBA2+)
# C6: Naive / Memory CD4+ T cells (CD3D+, IL7R+, CCR7+)
# C7: B cells (MS4A1/CD20+, CD79A+, CD79B+, BANK1+)
# C8: Naive T cells (CD3D+, CCR7+, IL7R+, high ribosomal)
# C9: Plasmacytoid dendritic cells / pDC (TCF4+, CCDC50+, IRF7+, IRF8+, IL3RA+)
# C10: CD16+ Natural killer cells / Mature NK cells (NKG7+, PRF1+, GZMB+, FCGR3A+, KLRF1+)
# C11: Conventional dendritic cells 1 / cDC1 (CLEC9A+, C1orf54+, WDFY4+, CADM1+)
# C12: Hematopoietic stem and progenitor cells / HSPCs (CD34+, SOX4+, CDK6+)
# C13: Classical monocytes (CD14+, S100A9+, S100A8+, VCAN+, LYZ+)
# C14: Effector CD8+ T cells (CD3D+, CD8A+, NKG7+, CCL5+, GZMA+, GZMH+, CST7+)
# C15: Proliferating NK cells (MKI67+, TOP2A+, NKG7+, PRF1+, GZMB+)
# C16: Innate lymphoid cells / ILCs (KLRB1+, IL7R+, GATA3+, TNFRSF18+, DDIT4+)
# C17: CD56bright Natural killer cells / NK cells (NCAM1+, KLRC1+, XCL1+, XCL2+, KLRD1+, GNLY+)
# C18: Plasma cells (MZB1+, JCHAIN+, SDC1-, TNFRSF17+)

cluster_to_celltype = {
    '0': 'Proliferating T cells',
    '1': 'Non-classical monocytes',
    '2': 'Platelets',
    '3': 'cDC2',
    '4': 'Memory CD8+ T cells',
    '5': 'Erythrocytes',
    '6': 'CD4+ T cells',
    '7': 'B cells',
    '8': 'Naive T cells',
    '9': 'Plasmacytoid dendritic cells',
    '10': 'NK cells',
    '11': 'cDC1',
    '12': 'Hematopoietic stem and progenitor cells',
    '13': 'Classical monocytes',
    '14': 'Effector CD8+ T cells',
    '15': 'Proliferating NK cells',
    '16': 'Innate lymphoid cells',
    '17': 'NK cells',
    '18': 'Plasma cells'
}

raw_adata = sc.read_h5ad('raw_counts.h5ad')
df_labels = pd.DataFrame({
    'barcode': raw_adata.obs_names,
    'cell_type': adata.obs['leiden'].map(cluster_to_celltype).values
})

df_labels.to_csv('labels.csv', index=False)
print("labels.csv written successfully! Shape:", df_labels.shape)
print("Distribution:\n", df_labels['cell_type'].value_counts())