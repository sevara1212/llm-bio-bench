import scanpy as sc

adata = sc.read_h5ad('raw_counts.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

genes_to_check = [
    'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B',  # T cells
    'IL7R', 'CCR7', 'S100A4', 'FOXP3', 'CD25', # CD4 T, Naive/Memory, Treg
    'GZMK', 'GZMB', 'GZMA', 'PRF1', 'NKG7', 'GNLY', 'NCAM1', 'FCGR3A', # Cytotoxic, NK
    'MS4A1', 'CD19', 'CD79A', 'CD79B', 'MZB1', 'SDC1', # B cells, Plasma
    'CD14', 'FCGR3A', 'LYZ', 'S100A9', 'S100A8', 'MS4A7', # Monocytes (classical vs non-classical)
    'CST3', 'HLA-DRA', 'CLEC9A', 'CD1C', 'IL3RA', 'CLEC4C', 'LILRA4', 'TCF4', # Dendritic cells (cDC1, cDC2, pDC)
    'PPBP', 'PF4', # Megakaryocytes / Platelets
    'HBB', 'HBA1', # Erythrocytes
    'MKI67', 'TOP2A', 'STMN1', # Proliferating
    'CD34', 'SOX4', 'PRSS57' # Progenitor / HSPC
]

# Check presence
genes_present = [g for g in genes_to_check if g in adata.var_names]
print(f"Present genes ({len(genes_present)}):", genes_present)