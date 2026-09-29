import scanpy as sc, pandas as pd, numpy as np
ad=sc.read_h5ad('processed.h5ad'); k='leiden_0.4'
labels={
'0':'Naive CD4 T cell','1':'B cell','2':'Cytotoxic CD8 T cell','3':'Memory CD4 T cell',
'4':'Proliferating NK cell','5':'CD161+ memory T cell','6':'CD14 monocyte','7':'FCGR3A monocyte',
'8':'NK cell','9':'Plasmacytoid dendritic cell','10':'Plasma cell','11':'Hematopoietic progenitor cell',
'12':'Platelet','13':'CD1C dendritic cell','14':'Cytotoxic NK cell','15':'CLEC9A dendritic cell',
'16':'GZMK+ CD8 T cell','17':'Gamma delta T cell','18':'Erythroid cell','19':'Plasmacytoid dendritic cell'}
ad.obs['cell_type']=ad.obs[k].map(labels).astype(str)
print(pd.crosstab(ad.obs[k],ad.obs.cell_type))
# All input barcodes, including three QC-failing barcodes explicitly designated low quality
orig=sc.read_h5ad('raw_counts.h5ad')
out=pd.DataFrame({'barcode':orig.obs_names.astype(str),'cell_type':'Low-quality cell'})
out.loc[out.barcode.isin(ad.obs_names), 'cell_type']=out.loc[out.barcode.isin(ad.obs_names),'barcode'].map(ad.obs.cell_type)
out.to_csv('labels.csv',index=False)
print('\noutput:',out.shape,'unique barcodes',out.barcode.nunique(),'missing labels',out.cell_type.isna().sum())
print(out.cell_type.value_counts())
print(out.head())
# persist annotation
ad.write_h5ad('processed_annotated.h5ad')