# Let's inspect canonical markers across clusters in leiden_0.3 and leiden_0.5
canonical_markers = [
    'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', # T cells
    'NCAM1', 'KLRD1', 'NKG7', 'GNLY', # NK cells
    'MS4A1', 'CD19', 'CD79A', # B cells
    'MZB1', 'SDC1', 'JCHAIN', # Plasma cells
    'CD14', 'LYZ', 'FCGR3A', 'MS4A7', # Monocytes (CD14+, CD16+)
    'PPBP', 'PF4', # Platelets / Megakaryocytes
    'IL3RA', 'CLEC4C', 'TCF4', 'LILRA4', # pDC
    'CLEC9A', 'CD1C', 'HLA-DRA', # cDC
    'HBA1', 'HBB', # Erythrocytes / RBCs
    'MKI67', 'TOP2A', 'STMN1', # Proliferating cells
    'CD34', 'SOX4' # Progenitors / HSC
]
present_markers = [m for m in canonical_markers if m in adata_proc.var_names]
print("Present canonical markers:", present_markers)

# Dotplot or matrixplot of these markers for leiden_0.3
import matplotlib.pyplot as plt
sc.pl.dotplot(adata_proc, var_names=present_markers, groupby='leiden_0.3', save='_leiden03.png')
sc.pl.dotplot(adata_proc, var_names=present_markers, groupby='leiden_0.5', save='_leiden05.png')
print("Plots saved.")