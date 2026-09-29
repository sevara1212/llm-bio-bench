import scanpy as sc
import pandas as pd

ad = sc.read_h5ad('processed.h5ad')

marker_sets = {
    'T_general': ['CD3D','CD3E','CD3G','TRAC'],
    'CD4_CD8': ['CD4','CD8A','CD8B'],
    'Tnaive_mem': ['IL7R','CCR7','SELL','TCF7'],
    'NK': ['NKG7','GNLY','KLRD1','KLRB1','GZMK','GZMB','FCGR3A'],
    'Mono': ['CD14','LYZ','S100A8','S100A9','FCN1','VCAN'],
    'B': ['MS4A1','CD79A','CD79B','CD19','IGHD','TCL1A','BANK1'],
    'Plasma': ['MZB1','JCHAIN','TNFRSF17'],
    'Platelet': ['PPBP','PF4','GNG11'],
    'RBC': ['HBB','HBA1'],
    'DC': ['LILRA4','IRF8','TCF4','CLEC9A','FCER1A','CD1C'],
    'Cycling': ['STMN1','MKI67','PCNA'],
    'Treg': ['FOXP3','IL2RA','CTLA4'],
}

for name, genes in marker_sets.items():
    genes = [g for g in genes if g in ad.raw.var_names]
    df = sc.get.obs_df(ad, keys=genes+['leiden'], use_raw=True)
    means = df.groupby('leiden').mean()
    print(f"--- {name} ---")
    print(means.round(2))
    print()