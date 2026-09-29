import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_final.h5ad')

markers_dict = {
    'T cell (general)': ['CD3D','CD3E','CD3G','TRAC','TRBC1'],
    'CD4 T': ['IL7R','CD4','CCR7','LEF1','SELL'],
    'CD8 T': ['CD8A','CD8B','GZMK','GZMH'],
    'Treg': ['FOXP3','IL2RA','TNFRSF18','CTLA4'],
    'NK': ['GNLY','NKG7','KLRD1','NCAM1','FCGR3A'],
    'B cell': ['MS4A1','CD79A','CD79B','CD19'],
    'Naive B': ['IGHD','TCL1A'],
    'Plasma': ['MZB1','JCHAIN','XBP1'],
    'Monocyte': ['CD14','LYZ','FCN1','CTSS','S100A8'],
    'cDC': ['FCER1A','CLEC10A','CD1C'],
    'cDC1': ['CLEC9A','BATF3','XCR1'],
    'pDC': ['LILRA4','IL3RA','CLEC4C','TCF4'],
    'Platelet': ['PPBP','PF4','TUBB1'],
    'HSPC': ['CD34','PRSS57','SOX4'],
    'Erythroid': ['HBB','HBA1','GYPA'],
    'Proliferating': ['MKI67','TOP2A','STMN1']
}

# use raw log-normalized data
adata_raw = adata.raw.to_adata()
adata_raw.obs['leiden'] = adata.obs['leiden'].values

for celltype, genes in markers_dict.items():
    genes_present = [g for g in genes if g in adata_raw.var_names]
    if not genes_present:
        continue
    df = pd.DataFrame(adata_raw[:, genes_present].X.toarray(), columns=genes_present, index=adata_raw.obs_names)
    df['leiden'] = adata_raw.obs['leiden'].values
    means = df.groupby('leiden').mean()
    print(f"\n=== {celltype} ===")
    print(means.round(2))