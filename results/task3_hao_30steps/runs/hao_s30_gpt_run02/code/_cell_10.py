import scanpy as sc, pandas as pd, numpy as np
ad=sc.read_h5ad('processed_markers.h5ad'); key='leiden_0.6'
labels={
'0':'Naive CD4 T cell','1':'Memory CD4 T cell','2':'Memory B cell','3':'Cytotoxic T cell',
'4':'Cycling lymphocyte','5':'CD16 NK cell','6':'Naive CD8 T cell','7':'CD161+ memory T cell',
'8':'Conventional dendritic cell','9':'CD56bright NK cell','10':'CD16 monocyte','11':'CD14 monocyte',
'12':'Hematopoietic progenitor cell','13':'Plasmacytoid dendritic cell','14':'Plasma cell','15':'Naive B cell',
'16':'Platelet','17':'GZMK+ memory T cell','18':'CD1C- conventional dendritic cell','19':'Gamma delta T cell',
'20':'Erythroid cell','21':'Plasmacytoid dendritic cell'}
# Typo check marker refined cluster 18 CLEC9A is cDC1, amend
labels['18']='CD141+ conventional dendritic cell'
# all input barcodes, mark the only QC exclusions transparently
raw=sc.read_h5ad('raw_counts.h5ad')
out=pd.DataFrame({'barcode':raw.obs_names.astype(str)})
assigned=pd.Series(ad.obs[key].astype(str).map(labels).values,index=ad.obs_names)
out['cell_type']=out.barcode.map(assigned).fillna('Low-quality cell')
out.to_csv('labels.csv',index=False)
# write non-housekeeping ranked markers for review
badprefix=('RPL','RPS','MT-','MALAT1','B2M','TMSB','EEF','GAPDH','ACTB','HBA','HBB')
rows=[]
for cl in ad.obs[key].cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=cl)
 d=d[~d['names'].str.startswith(badprefix)].head(15)
 for rank,(_,r) in enumerate(d.iterrows(),1):
  rows.append({'cluster':cl,'cell_type':labels[str(cl)],'n_cells':int((ad.obs[key]==cl).sum()),'rank':rank,'gene':r['names'],'score':r['scores'],'logfoldchange':r['logfoldchanges'],'pval_adj':r['pvals_adj']})
pd.DataFrame(rows).to_csv('markers.csv',index=False)
print(out.shape, out.cell_type.value_counts().to_string())
print('missing',out.isna().sum().to_dict(),'duplicates',out.barcode.duplicated().sum())
print('marker rows',len(rows)); print(pd.read_csv('labels.csv').head())