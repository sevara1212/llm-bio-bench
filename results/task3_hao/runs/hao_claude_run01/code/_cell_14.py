import scanpy as sc
import numpy as np
import pandas as pd

ad = sc.read_h5ad('processed.h5ad')

# Comprehensive marker sets for score_genes
signatures = {
    'CD14+ Monocyte': ['CD14','LYZ','S100A8','S100A9','FCN1','VCAN','S100A12','MNDA','CSF3R'],
    'CD16+ Monocyte': ['FCGR3A','MS4A7','LST1','AIF1','CFD','SERPINA1','COTL1'],
    'CD4 T Naive': ['CD3D','CD3E','TRAC','CD4','CCR7','SELL','TCF7','LEF1'],
    'CD4 T Memory': ['CD3D','CD3E','TRAC','CD4','IL7R','KLRB1'],
    'CD8 T Naive': ['CD3D','CD3E','TRAC','CD8A','CD8B','CCR7','SELL','TCF7'],
    'CD8 T Effector/Memory': ['CD3D','CD3E','TRAC','CD8A','CD8B','GZMK','GZMA','CD27'],
    'NK': ['NKG7','GNLY','KLRD1','KLRF1','FCGR3A','GZMB','PRF1','XCL1','XCL2'],
    'B naive': ['MS4A1','CD79A','CD79B','CD19','IGHD','TCL1A'],
    'B memory': ['MS4A1','CD79A','CD79B','BANK1','CD27'],
    'Plasma cell': ['MZB1','JCHAIN','TNFRSF17','SEC11C','DERL3'],
    'Platelet': ['PPBP','PF4','GNG11','TUBB1','NRGN'],
    'RBC': ['HBB','HBA1','HBA2','ALAS2'],
    'pDC': ['LILRA4','IRF8','TCF4','IL3RA','BCL11A','SERPINF1','PLD4'],
    'cDC': ['CD1C','FCER1A','CLEC9A','BATF3','XCR1','CD74','HLA-DRA'],
    'Cycling/Progenitor': ['STMN1','MKI67','PCNA','CD34','PRSS57','SOX4'],
}

for name, genes in signatures.items():
    genes = [g for g in genes if g in ad.raw.var_names]
    sc.tl.score_genes(ad, genes, score_name=f'score_{name}', use_raw=True)

score_cols = [f'score_{n}' for n in signatures]
df_scores = ad.obs[score_cols + ['leiden']].groupby('leiden').mean()
df_scores.columns = list(signatures.keys())
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)
print(df_scores.round(3))

# best type per cluster
best = df_scores.idxmax(axis=1)
print(best)