import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('full_with_clusters.h5ad')

marker_genes = {
    'T cell (general)': ['CD3D','CD3E','CD3G','TRAC','TRBC2'],
    'CD4 T': ['CD4','IL7R'],
    'CD8 T': ['CD8A','CD8B'],
    'Naive T': ['CCR7','SELL','TCF7','LEF1'],
    'Treg': ['FOXP3','IL2RA'],
    'NK': ['GNLY','NKG7','KLRD1','KLRB1','NCAM1','FCGR3A','KLRF1'],
    'B cell': ['MS4A1','CD79A','CD79B','CD19','BANK1'],
    'Naive B': ['TCL1A','IGHD','IGHM'],
    'Plasma': ['MZB1','JCHAIN','TNFRSF17','SDC1'],
    'Monocyte CD14': ['CD14','LYZ','S100A8','S100A9','VCAN','FCN1'],
    'Monocyte CD16': ['FCGR3A','MS4A7','LST1','CDKN1C'],
    'cDC': ['CST3','FCER1A','CLEC9A','CD1C','ITGAX'],
    'pDC': ['LILRA4','IRF7','TCF4','IL3RA','CLEC4C'],
    'Platelet': ['PPBP','PF4','GNG11','ITGA2B'],
    'Erythrocyte': ['HBB','HBA1','HBA2','ALAS2'],
    'Proliferating': ['MKI67','STMN1','HMGB2','TOP2A'],
    'Progenitor/HSPC': ['CD34','PRSS57','SOX4'],
    'Cytotoxic': ['GZMB','GZMA','GZMK','GZMH','PRF1','NKG7','CST7'],
}

df = pd.DataFrame(index=sorted(adata.obs.leiden.unique(), key=int))
for cell_type, genes in marker_genes.items():
    genes_present = [g for g in genes if g in adata.var_names]
    if not genes_present:
        continue
    expr = adata[:, genes_present].X
    if not isinstance(expr, np.ndarray):
        expr = expr.toarray()
    mean_expr = expr.mean(axis=1)
    adata.obs[f'score_{cell_type}'] = mean_expr
    df[cell_type] = adata.obs.groupby('leiden')[f'score_{cell_type}'].mean()

pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 20)
print(df.round(2))