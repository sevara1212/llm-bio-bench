# Let's check expression of canonical markers across these clusters
canonical_markers = [
    'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B',  # T cells
    'IL7R', 'CCR7', 'S100A4',                # CD4 naive/memory
    'GZMA', 'GZMB', 'GZMK', 'NKG7', 'PRF1',  # Cytotoxic / NK
    'NCAM1', 'FCGR3A',                       # CD56, CD16
    'MS4A1', 'CD19', 'CD79A',                # B cells
    'MZB1', 'SDC1', 'CD38',                  # Plasma cells
    'CD14', 'LYZ', 'S100A9',                 # Classical monocytes
    'FCGR3A', 'MS4A7',                       # Non-classical monocytes
    'FCER1A', 'CST3', 'CLEC9A',              # cDCs
    'IL3RA', 'CLEC4C', 'TCF4',               # pDCs
    'PPBP', 'PF4',                           # Megakaryocytes / Platelets
    'HBB', 'HBA1',                           # Erythrocytes
    'MKI67', 'TOP2A', 'STMN1',               # Proliferating
    'CD34'                                   # Progenitors / HSCs
]
available_markers = [m for m in canonical_markers if m in adata.raw.var_names]

# Let's compute mean expression per cluster for these markers
mean_expr = pd.DataFrame(index=adata.obs['leiden_0.5'].cat.categories)
for m in available_markers:
    mean_expr[m] = [adata.raw[adata.obs['leiden_0.5'] == cl, m].X.mean() for cl in mean_expr.index]

print(mean_expr[['CD3D', 'CD4', 'CD8A', 'MS4A1', 'CD14', 'FCGR3A', 'PPBP', 'HBB', 'MZB1', 'TCF4', 'CLEC9A', 'STMN1']])