import scanpy as sc
import pandas as pd
import numpy as np

full = sc.read_h5ad('full_with_markers.h5ad')

markers = {
    'CD3D':'T cell','CD3E':'T cell','CD3G':'T cell','TRAC':'T cell',
    'CD4':'CD4 T','IL7R':'CD4 T','CCR7':'Naive T',
    'CD8A':'CD8 T','CD8B':'CD8 T',
    'GZMK':'CD8 memory','GZMB':'Cytotoxic','GZMH':'Cytotoxic','NKG7':'NK/CTL','GNLY':'NK','KLRD1':'NK','KLRF1':'NK','NCAM1':'NK',
    'MS4A1':'B cell','CD79A':'B cell','CD79B':'B cell','IGHD':'Naive B','IGHM':'B','TCL1A':'Naive B','BANK1':'B',
    'MZB1':'Plasma','JCHAIN':'Plasma','TNFRSF17':'Plasma',
    'CD14':'Monocyte','LYZ':'Monocyte','S100A8':'Monocyte','S100A9':'Monocyte','FCN1':'Monocyte',
    'FCGR3A':'Mono(CD16)','MS4A7':'Mono(CD16)',
    'CST3':'DC','CD74':'DC/APC','FCER1G':'Myeloid',
    'CLEC9A':'cDC1','BATF3':'cDC1',
    'LILRA4':'pDC','IL3RA':'pDC','IRF7':'pDC',
    'PPBP':'Platelet','PF4':'Platelet','TUBB1':'Platelet',
    'HBB':'RBC','HBA1':'RBC','ALAS2':'RBC',
    'STMN1':'Proliferating','MKI67':'Proliferating',
    'PRSS57':'HSPC','CD34':'HSPC',
    'FOXP3':'Treg','TNFRSF18':'Treg/Act',
}
genes_present = [g for g in markers if g in full.var_names]
df = sc.get.obs_df(full, keys=genes_present+['leiden'])
mean_expr = df.groupby('leiden').mean()
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 50)
print(mean_expr.round(2).T)