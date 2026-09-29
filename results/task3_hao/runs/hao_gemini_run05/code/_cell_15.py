# Let's inspect CD3D, CD4, CD8A, CD8B, CCR7, LEF1 for clusters 0, 3, 6, 12, 15, 17
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('processed.h5ad')

check_genes = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'LEF1', 'SELL', 'IL7R', 'FOXP3', 'IKZF2', 'TRDC', 'TRGC1']
sub_adata = adata[:, [g for g in check_genes if g in adata.var_names]]

for c in ['0', '2', '3', '6', '8', '12', '14', '15', '17']:
    c_data = sub_adata[sub_adata.obs['leiden'] == c]
    means = np.asarray(c_data.X.mean(axis=0)).flatten()
    print(f"Cluster {c}: " + ", ".join([f"{g}:{val:.2f}" for g, val in zip(sub_adata.var_names, means)]))

# Cluster annotations:
# Cluster 0: CD3D:0.77, CD4:0.35, CD8A:0.04, CCR7:0.43, LEF1:0.38, SELL:0.67, IL7R:0.62 -> Naive CD4+ T cell
# Cluster 1: MS4A1, CD79A, BANK1, CD79B -> B cell
# Cluster 2: CD3D:0.75, CD8A:0.87, CD8B:0.62, NKG7, CCL5, GZMH, GZMA -> Effector CD8+ T cell / Cytotoxic CD8+ T cell
# Cluster 3: CD3D:0.84, CD4:0.35, CD8A:0.04, CCR7:0.12, IL7R:0.91, S100A4/IL32 -> Memory CD4+ T cell
# Cluster 4: STMN1, DEK, HMGB2, TYMS, PCNA, PCLAF -> Proliferating cell
# Cluster 5: VCAN, CD14, FCN1, S100A8, S100A9, LYZ -> CD14+ Monocyte
# Cluster 6: CD3D, KLRB1 (CD161), DUSP2, GZMK, RORA, NCR3 -> CD8+ / MAIT / Memory T cell
# Cluster 7: LST1, AIF1, CFD, MS4A7, FCGR3A -> FCGR3A+ / CD16+ Monocyte
# Cluster 8: KLRD1, GNLY, XCL2, KLRF1, KLRC1, IL2RB, CD3D:0.05 -> NK cell (CD56bright / NK cell)
# Cluster 9: HLA-DRA, CST3, FCER1A, HLA-DPA1 -> Conventional Dendritic Cell (cDC2)
# Cluster 10: CCDC50, TCF4, PLD4, IRF7, IRF8, IL3RA (CD123), LILRA4 -> Plasmacytoid Dendritic Cell (pDC)
# Cluster 11: MZB1, HSP90B1, JCHAIN, TNFRSF17 (BCMA), SDC1 -> Plasma cell
# Cluster 12: CD3D:0.92, CD4:0.41, CDK6, SOX4, RPS24 -> Naive CD4+ T cell (or T cell)
# Cluster 13: PF4, PPBP, CAVIN2, TUBB1 -> Platelet / Megakaryocyte
# Cluster 14: PRF1, NKG7, CST7, KLRF1, KLRD1, GZMB, FCGR3A, CD3D:0.06 -> NK cell (CD56dim NK cell)
# Cluster 15: CD3D:0.96, CD8A:0.43, GZMK, CD27, TCF7 -> Central memory / GZMK+ CD8+ T cell
# Cluster 16: CD74, HLA-DPA1, CLEC9A, WDFY4, BATF3 -> Conventional Dendritic Cell (cDC1)
# Cluster 17: TNFRSF18 (GITR), TNFRSF4 (OX40), FOXP3, T cell -> Regulatory T cell (Treg)
# Cluster 18: LILRA4, TCF4, PLD4, SAMHD1, IRF8 -> Plasmacytoid Dendritic Cell (pDC)

# Create mapping dictionary and write labels.csv
cluster_labels = {
    '0': 'Naive CD4+ T cell',
    '1': 'B cell',
    '2': 'Cytotoxic CD8+ T cell',
    '3': 'Memory CD4+ T cell',
    '4': 'Proliferating cell',
    '5': 'CD14+ Monocyte',
    '6': 'Effector memory T cell',
    '7': 'CD16+ Monocyte',
    '8': 'NK cell',
    '9': 'Conventional dendritic cell',
    '10': 'Plasmacytoid dendritic cell',
    '11': 'Plasma cell',
    '12': 'Naive CD4+ T cell',
    '13': 'Platelet',
    '14': 'NK cell',
    '15': 'Memory CD8+ T cell',
    '16': 'Conventional dendritic cell',
    '17': 'Regulatory T cell',
    '18': 'Plasmacytoid dendritic cell'
}

raw_adata = sc.read_h5ad('raw_counts.h5ad')
cell_type_series = adata.obs['leiden'].map(cluster_labels)
out_df = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': cell_type_series.values
})
out_df.to_csv('labels.csv', index=False)
print(f"Saved {len(out_df)} labels to labels.csv")
print("Cell type distribution:")
print(out_df['cell_type'].value_counts())
print("\nCheck missing barcodes from raw_counts:")
print("Matches raw_counts exactly:", (out_df['barcode'] == raw_adata.obs_names).all())