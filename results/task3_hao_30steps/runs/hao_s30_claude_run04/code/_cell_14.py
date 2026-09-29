import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_full_annotated.h5ad')

markers = {
    'CD3D':'T cell', 'CD3E':'T cell', 'CD4':'CD4 T', 'CD8A':'CD8 T', 'CD8B':'CD8 T',
    'IL7R':'CD4 T', 'CCR7':'Naive T', 'SELL':'Naive T', 'TCF7':'Naive T',
    'GZMK':'Memory CD8 T', 'GZMB':'Cytotoxic', 'GZMH':'Cytotoxic', 'NKG7':'NK/CTL', 'GNLY':'NK', 'KLRD1':'NK',
    'FCGR3A':'NK/Mono', 'NCAM1':'NK',
    'MS4A1':'B cell', 'CD79A':'B cell', 'CD79B':'B cell', 'CD19':'B cell', 'IGHD':'Naive B', 'TCL1A':'Naive B',
    'MZB1':'Plasma', 'JCHAIN':'Plasma', 'TNFRSF17':'Plasma',
    'CD14':'Mono', 'LYZ':'Mono', 'S100A8':'Mono', 'FCN1':'Mono', 'VCAN':'Mono',
    'FCGR3A':'Mono CD16', 'MS4A7':'Mono CD16', 'LST1':'Mono CD16',
    'CST3':'DC', 'FCER1A':'DC', 'CLEC9A':'cDC1', 'CD1C':'cDC2',
    'IL3RA':'pDC', 'IRF7':'pDC', 'IRF8':'pDC',
    'PPBP':'Platelet', 'PF4':'Platelet', 'GNG11':'Platelet',
    'HBB':'RBC', 'HBA1':'RBC', 'HBA2':'RBC',
    'FOXP3':'Treg', 'IL2RA':'Treg', 'TNFRSF18':'Treg', 'TNFRSF4':'Treg',
    'STMN1':'Proliferating', 'MKI67':'Proliferating',
    'PRSS57':'Progenitor', 'CD34':'Progenitor', 'SOX4':'Progenitor'
}
genes_present = [g for g in markers if g in adata.var_names]
print(len(genes_present), 'of', len(markers))
missing = [g for g in markers if g not in adata.var_names]
print("missing:", missing)