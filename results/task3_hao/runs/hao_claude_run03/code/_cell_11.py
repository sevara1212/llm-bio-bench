import scanpy as sc
import pandas as pd

full = sc.read_h5ad('full_with_markers.h5ad')

markers = ['CD3D','CD3E','TRAC','IL7R','CCR7','CD4','CD8A','CD8B','GZMK','GZMH','GZMB','NKG7','GNLY','KLRD1','KLRF1',
           'MS4A1','CD79A','IGHD','TCL1A','BANK1','MZB1','JCHAIN','TNFRSF17',
           'CD14','LYZ','S100A8','S100A9','FCN1','FCGR3A','MS4A7','CST3','CD74','FCER1G',
           'CLEC9A','BATF3','LILRA4','IL3RA','IRF7','PPBP','PF4','TUBB1','HBB','HBA1','ALAS2',
           'STMN1','MKI67','PRSS57','CD34','FOXP3','TNFRSF18']
genes_present = [g for g in markers if g in full.var_names]
df = sc.get.obs_df(full, keys=genes_present+['leiden'])
mean_expr = df.groupby('leiden').mean()

for cl in mean_expr.index:
    row = mean_expr.loc[cl].sort_values(ascending=False)
    print(cl, row.head(8).round(2).to_dict())