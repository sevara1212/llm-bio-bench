import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('adata_processed.h5ad')

markers = {
 'CD3D':'T', 'CD3E':'T','CD3G':'T','TRAC':'T','IL7R':'T_CD4',
 'CD4':'T_CD4','CD8A':'T_CD8','CD8B':'T_CD8','CCR7':'Naive_T',
 'GZMK':'CD8_Tem','GZMB':'CD8_Tem','NKG7':'NK','GNLY':'NK','KLRD1':'NK','KLRF1':'NK','NCAM1':'NK',
 'MS4A1':'B','CD79A':'B','CD79B':'B','CD19':'B','BANK1':'B',
 'CD14':'Mono','LYZ':'Mono','S100A8':'Mono_classical','S100A9':'Mono_classical','FCN1':'Mono_classical',
 'FCGR3A':'Mono_nonclassical','MS4A7':'Mono_nonclassical',
 'CST3':'DC','FCER1A':'cDC','CLEC10A':'cDC','CLEC9A':'cDC1',
 'LILRA4':'pDC','IL3RA':'pDC','CLEC4C':'pDC',
 'PPBP':'Platelet','PF4':'Platelet','GP9':'Platelet',
 'HBB':'Erythrocyte','HBA1':'Erythrocyte','HBA2':'Erythrocyte',
 'MZB1':'Plasma','JCHAIN':'Plasma','IGHG1':'Plasma',
 'MKI67':'Proliferating','STMN1':'Proliferating','TYMS':'Proliferating',
 'FOXP3':'Treg','IL2RA':'Treg'
}

genes_present = [g for g in markers if g in adata.raw.var_names]
print(genes_present)

import numpy as np
expr = adata.raw[:, genes_present].X
import scipy.sparse as sp
if sp.issparse(expr):
    expr = expr.toarray()
df = pd.DataFrame(expr, columns=genes_present, index=adata.obs_names)
df['leiden'] = adata.obs['leiden'].values
mean_expr = df.groupby('leiden').mean()
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)
print(mean_expr.T)