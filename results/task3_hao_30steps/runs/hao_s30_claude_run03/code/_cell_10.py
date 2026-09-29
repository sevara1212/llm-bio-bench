import scanpy as sc
import pandas as pd
import numpy as np

ad = sc.read_h5ad('clustered.h5ad')

markers = {
    'CD3D':'T', 'CD3E':'T', 'CD4':'CD4T', 'CD8A':'CD8T', 'CD8B':'CD8T',
    'IL7R':'CD4T_naive/memory', 'CCR7':'Naive', 'SELL':'Naive',
    'GZMK':'Tem', 'GZMB':'Cytotoxic', 'GZMH':'Cytotoxic', 'NKG7':'NK/CTL', 'GNLY':'NK',
    'KLRD1':'NK', 'KLRF1':'NK', 'FCGR3A':'NK/Mono',
    'MS4A1':'B', 'CD79A':'B', 'CD79B':'B', 'IGHD':'Naive B', 'TCL1A':'Naive B',
    'MZB1':'Plasma', 'JCHAIN':'Plasma', 'TNFRSF17':'Plasma',
    'CD14':'CD14Mono', 'LYZ':'Mono', 'S100A8':'CD14Mono', 'S100A9':'CD14Mono', 'FCN1':'CD14Mono',
    'MS4A7':'CD16Mono', 'LST1':'CD16Mono', 'AIF1':'Mono/DC',
    'CST3':'DC/Mono', 'FCER1A':'cDC2', 'CD1C':'cDC2', 'CLEC9A':'cDC1', 'BATF3':'cDC1',
    'LILRA4':'pDC', 'IL3RA':'pDC', 'IRF7':'pDC', 'TCF4':'pDC',
    'PPBP':'Platelet', 'PF4':'Platelet', 'NRGN':'Platelet',
    'HBB':'Erythrocyte', 'HBA1':'Erythrocyte',
    'STMN1':'Proliferating', 'MKI67':'Proliferating', 'PRSS57':'HSPC', 'CD34':'HSPC'
}

df = sc.get.obs_df(ad, keys=list(markers.keys())+['leiden_05'], use_raw=True)
mean_expr = df.groupby('leiden_05').mean()
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 50)
print(mean_expr.round(2))