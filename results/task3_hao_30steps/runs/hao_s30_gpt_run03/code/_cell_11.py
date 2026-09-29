import scanpy as sc, numpy as np, pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from scipy.sparse import issparse
ad=sc.read_h5ad('processed.h5ad'); raw=sc.read_h5ad('raw_counts.h5ad'); key='leiden_0.6'
labels={'0':'Naive CD4 T cell','1':'Memory CD4 T cell','2':'Effector memory CD8 T cell','3':'Memory B cell','4':'Cycling NK cell','5':'CD161+ memory T cell','6':'CD16+ NK cell','7':'CD1c+ conventional dendritic cell','8':'CD56bright NK cell','9':'Naive B cell','10':'CD16+ monocyte','11':'Plasmacytoid dendritic cell','12':'CD14+ monocyte','13':'Mast cell','14':'Plasma cell','15':'Platelet','16':'CLEC9A+ dendritic cell','17':'Naive CD8 T cell','18':'Activated T cell','19':'Plasmacytoid dendritic cell'}
ad.obs['cell_type']=ad.obs[key].map(labels).astype(str)
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=100,key_added='markers')
sc.get.rank_genes_groups_df(ad,group=None,key='markers').to_csv('cluster_markers.csv',index=False)
X=raw[:,ad.var_names].X.copy().astype(float).tocsr(); lib=np.asarray(X.sum(axis=1)).ravel(); lib[lib==0]=1
X=X.multiply((1e4/lib)[:,None]).tocsr(); X.data=np.log1p(X.data)
idx=np.flatnonzero(ad.var.highly_variable.values); Z=X[:,idx].toarray();Z=(Z-ad.var['mean'].values[idx])/ad.var['std'].values[idx];Z=np.clip(Z,-10,10)
P=Z@ad.varm['PCs'][idx,:]
knn=KNeighborsClassifier(n_neighbors=25,weights='distance',n_jobs=1).fit(ad.obsm['X_pca'][:,:40],ad.obs.cell_type)
pred=knn.predict(P[:,:40])
gene=np.asarray(raw.var_names); hb=np.array([g.startswith('HBA') or g.startswith('HBB') or g in ('HBD','HBM','HBQ1') for g in gene]); hbfrac=np.asarray(raw[:,hb].X.sum(axis=1)).ravel()/lib
pred[hbfrac>=.20]='Erythrocyte'
out=pd.DataFrame({'barcode':raw.obs_names.astype(str),'cell_type':pred});lookup=ad.obs.cell_type.to_dict();out['cell_type']=[lookup.get(b,ct) for b,ct in zip(out.barcode,out.cell_type)];out.to_csv('labels.csv',index=False)
ad.write_h5ad('processed_annotated.h5ad')
print('labels rows/unique',out.shape);print(out.cell_type.value_counts().to_string()); print('missing',out.isna().sum().to_dict(),'retained',sum(out.barcode.isin(ad.obs_names)),'excluded',len(raw)-len(ad),'RBC',sum(hbfrac>=.2))