import scanpy as sc
import numpy as np
import pandas as pd

ad = sc.read_h5ad('processed.h5ad')

markers = ['CD3D','CD3E','CD3G','TRAC','CD4','CD8A','CD8B','IL7R','CCR7','SELL','TCF7',
           'NKG7','GNLY','KLRD1','KLRB1','GZMK','GZMB','FCGR3A',
           'CD14','LYZ','S100A8','S100A9','FCN1','VCAN',
           'MS4A1','CD79A','CD79B','CD19','IGHD','TCL1A','BANK1',
           'MZB1','JCHAIN','TNFRSF17',
           'PPBP','PF4','GNG11',
           'HBB','HBA1',
           'LILRA4','IRF8','TCF4','CLEC9A','FCER1A','CD1C',
           'STMN1','MKI67','PCNA',
           'FOXP3','IL2RA','CTLA4']

markers = [m for m in markers if m in ad.raw.var_names]
df = sc.get.obs_df(ad, keys=markers+['leiden'], use_raw=True)
means = df.groupby('leiden').mean()
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
print(means.round(2))