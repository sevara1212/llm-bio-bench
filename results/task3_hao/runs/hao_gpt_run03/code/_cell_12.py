import scanpy as sc, pandas as pd, numpy as np
from scipy.sparse import issparse
ad=sc.read_h5ad('processed.h5ad')
key='leiden_0.5'
labels={
'0':'Naive CD4 T cell','1':'Memory CD4 T cell','2':'B cell','3':'Cytotoxic T cell',
'4':'Cycling NK cell','5':'CD161+ memory T cell','6':'Conventional dendritic cell','7':'FCGR3A+ monocyte',
'8':'CD56bright NK cell','9':'Plasmacytoid dendritic cell','10':'Plasmablast','11':'Hematopoietic progenitor cell',
'12':'CD14+ monocyte','13':'Platelet','14':'CD56dim NK cell','15':'CD1C- conventional dendritic cell',
'16':'Activated CD4 T cell','17':'Memory T cell','18':'Erythrocyte','19':'Plasmacytoid dendritic cell'}
ad.obs['cell_type']=ad.obs[key].map(labels).astype(str)
print(pd.crosstab(ad.obs[key],ad.obs.cell_type))
# marker table export calculated rank stats, avoids h5 save
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=50)
rows=[]; r=ad.uns['rank_genes_groups']
for c in ad.obs[key].cat.categories:
 for rank,g in enumerate(r['names'][c][:50],1): rows.append([c,labels[c],rank,g,float(r['scores'][c][rank-1]),float(r['pvals_adj'][c][rank-1])])
pd.DataFrame(rows,columns=['cluster','cell_type','rank','gene','score','pval_adj']).to_csv('cluster_markers.csv',index=False)
# save compact annotations
ad.obs[[key,'cell_type','total_counts','n_genes_by_counts','pct_counts_mt','batch']].to_csv('cell_qc_cluster_labels.csv')
# label four QC-excluded cells using nearest centroid on library-size-normalised log counts, cosine dot product
raw=sc.read_h5ad('raw_counts.h5ad'); raw.var['mt']=raw.var_names.str.upper().str.startswith('MT-');sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True,log1p=False)
keep=(raw.obs.n_genes_by_counts>=500)&(raw.obs.pct_counts_mt<15)
print('kept',keep.sum(),'excluded',list(raw.obs_names[~keep]))
# Make labels indexed original
out=pd.Series(index=raw.obs_names,dtype=object); out.loc[ad.obs_names]=ad.obs.cell_type.values
# centroid based on marker gene log normalized values
sc.pp.normalize_total(raw,target_sum=1e4);sc.pp.log1p(raw)
# cosine on HVG available in processed var
features=ad.var_names.intersection(raw.var_names); fi=raw.var_names.get_indexer(features)
X=raw.X[:,fi]; X=X.toarray() if issparse(X) else X
X=X/(np.linalg.norm(X,axis=1,keepdims=True)+1e-12)
for bc in raw.obs_names[~keep]:
 i=raw.obs_names.get_loc(bc); bestlab=None; best=-9
 for lab in sorted(set(labels.values())):
  ids=raw.obs_names.get_indexer(ad.obs_names[ad.obs.cell_type==lab])
  cent=X[ids].mean(0); cent/=np.linalg.norm(cent)+1e-12
  score=X[i]@cent
  if score>best: best,bestlab=score,lab
 out.loc[bc]=bestlab; print(bc,bestlab,round(float(best),3))
pd.DataFrame({'barcode':raw.obs_names,'cell_type':out.loc[raw.obs_names].values}).to_csv('labels.csv',index=False)
print(pd.read_csv('labels.csv').head());print('rows',len(out),'missing',out.isna().sum());print(out.value_counts())