import scanpy as sc, pandas as pd, numpy as np
raw=sc.read_h5ad('raw_counts.h5ad')
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.6'
labels={
'0':'Naive CD4 T cell','1':'Cytotoxic T cell','2':'Memory CD4 T cell','3':'Memory B cell',
'4':'Cycling lymphocyte','5':'Conventional dendritic cell','6':'Naive CD8 T cell',
'7':'CD161+ memory T cell','8':'Non-classical monocyte','9':'Natural killer cell',
'10':'Naive B cell','11':'Plasmacytoid dendritic cell','12':'CD34+ hematopoietic progenitor',
'13':'Classical monocyte','14':'Plasma cell','15':'Natural killer cell','16':'Platelet',
'17':'Conventional dendritic cell 1','18':'Gamma delta T cell','19':'GZMK+ memory T cell','20':'Erythrocyte'}
assert set(ad.obs[key].cat.categories)==set(labels)
# markers table
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=100)
rows=[]
for c in ad.obs[key].cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=c)
 d=d[~d.names.str.startswith(('RPL','RPS','MT-'))].head(15).copy()
 d.insert(0,'cell_type',labels[c]);d.insert(0,'cluster',c)
 rows.append(d[['cluster','cell_type','names','scores','logfoldchanges','pvals_adj']])
pd.concat(rows).to_csv('marker_genes.csv',index=False)
# Map retained; explicitly retain all source barcodes, QC exclusions flagged
out=pd.DataFrame({'barcode':raw.obs_names.astype(str)})
maplab=pd.Series(ad.obs[key].astype(str).map(labels).values,index=ad.obs_names)
out['cell_type']=out.barcode.map(maplab).fillna('Low-quality cell')
out.to_csv('labels.csv',index=False)
print(out.cell_type.value_counts().to_string())
print('output',out.shape, 'unique barcodes',out.barcode.nunique(),'missing labels',out.cell_type.isna().sum())
print('\nCluster assignments:')
for c in ad.obs[key].cat.categories: print(c, labels[c], (ad.obs[key]==c).sum())
print('\nTop marker rows:');print(pd.read_csv('marker_genes.csv').groupby(['cluster','cell_type']).names.apply(lambda x:', '.join(x.head(5))).to_string())