import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('processed.h5ad')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

print("Top 10 marker genes per cluster:")
for g in groups:
    genes = [result['names'][g][i] for i in range(10)]
    print(f"Cluster {g}: {', '.join(genes)}")

# Also let's check expression of standard canonical PBMC markers across clusters:
markers = [
    'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', # T cells
    'IL7R', 'CCR7', 'S100A4', # CD4 memory/naive
    'GNLY', 'NKG7', 'NCAM1', 'FCGR3A', # NK / cytotoxic / CD16
    'MS4A1', 'CD19', 'CD79A', # B cells
    'CD14', 'LYZ', 'VCAN', 'S100A9', # CD14 Monocytes
    'FCGR3A', 'MS4A7', # FCGR3A+ / CD16+ Monocytes
    'FCER1A', 'CST3', 'CLEC10A', 'HLA-DPA1', # Dendritic cells (cDC)
    'LILRA4', 'IL3RA', 'CLEC4C', # pDC
    'PPBP', 'PF4', 'TUBB1', # Megakaryocytes / Platelets
    'HBB', 'HBA1', 'HBA2', # Erythrocytes
    'MZB1', 'IGKC', 'IGHG1', # Plasma cells
]
present_markers = [m for m in markers if m in adata.var_names]
print("\nPresent markers:", len(present_markers))