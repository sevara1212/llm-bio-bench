import pandas as pd,scanpy as sc, numpy as np
clusters=pd.read_csv('cluster_assignments_qc.csv',index_col=0).iloc[:,0].astype(str)
labels={
'0':'Cytotoxic T cell','1':'Naive CD4 T cell','2':'CD4 memory T cell','3':'Memory B cell','4':'Naive CD4 T cell',
'5':'Proliferating NK cell','6':'Conventional dendritic cell','7':'CD161+ memory T cell','8':'FCGR3A+ monocyte','9':'Natural killer cell',
'10':'Hematopoietic progenitor cell','11':'Naive B cell','12':'Plasmacytoid dendritic cell','13':'CD14+ monocyte','14':'Plasma cell',
'15':'Platelet','16':'Natural killer cell','17':'CLEC9A+ dendritic cell','18':'GZMK+ memory T cell','19':'Gamma delta T cell','20':'Erythroid cell'}
assert set(clusters.unique())==set(labels)
raw=sc.read_h5ad('raw_counts.h5ad')
out=pd.DataFrame({'barcode':raw.obs_names,'cell_type':'QC_filtered'})
out.loc[out.barcode.isin(clusters.index),'cell_type']=out.loc[out.barcode.isin(clusters.index),'barcode'].map(clusters.map(labels))
assert len(out)==raw.n_obs and out.barcode.is_unique and out.cell_type.notna().all()
out.to_csv('labels.csv',index=False)
summary=pd.DataFrame({'cluster':list(labels), 'cell_type':[labels[x] for x in labels], 'n_cells':[int((clusters==x).sum()) for x in labels]})
summary.to_csv('cluster_cell_type_summary.csv',index=False)
print(summary.to_string(index=False));print('\nlabels:',out.shape,out.cell_type.value_counts().to_dict());print(pd.read_csv('labels.csv').head());