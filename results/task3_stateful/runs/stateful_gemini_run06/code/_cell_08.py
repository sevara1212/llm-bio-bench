# Let's inspect known PBMC markers across the clusters
canonical_markers = {
    'CD4 T': ['IL7R', 'CD4', 'LDHB'],
    'CD8 T': ['CD8A', 'CD8B'],
    'Naive CD4 T': ['IL7R', 'CCR7', 'S100A4'],
    'Memory CD4 T': ['IL7R', 'S100A4'],
    'B cells': ['MS4A1', 'CD79A', 'CD79B'],
    'CD14+ Monocytes': ['CD14', 'LYZ', 'S100A9'],
    'FCGR3A+ Monocytes': ['FCGR3A', 'MS4A7', 'LST1', 'AIF1'],
    'NK cells': ['GNLY', 'NKG7'],
    'Dendritic cells': ['FCER1A', 'CST3', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DRA'],
    'Megakaryocytes / Platelets': ['PPBP', 'PF4', 'NRGN']
}

for group in marker_genes.columns:
    top_genes = marker_genes[group].tolist()
    print(f"Cluster {group}: {top_genes[:7]}")