import scanpy as sc, numpy as np, pandas as pd
from sklearn.neighbors import KNeighborsClassifier
ref=sc.read_h5ad('qc_processed.h5ad')
raw=sc.read_h5ad('raw_counts.h5ad')
# final cluster annotation, grounded in differential-marker review
mapping={'0':'Naive CD4 T cell','1':'B cell','2':'Cytotoxic T cell','3':'Memory CD4 T cell','4':'Cycling lymphocyte','5':'CD161+ memory T cell','6':'CD1c+ conventional dendritic cell','7':'FCGR3A+ monocyte','8':'Natural killer cell','9':'Plasmacytoid dendritic cell','10':'Hematopoietic stem/progenitor cell','11':'CD14+ monocyte','12':'Plasma cell','13':'Platelet','14':'Natural killer cell','15':'GZMK+ memory T cell','16':'CLEC9A+ conventional dendritic cell','17':'Gamma delta T cell','18':'Erythroid cell','19':'Plasmacytoid dendritic cell'}
ref.obs['cell_type']=ref.obs['leiden_0.5'].astype(str).map(mapping)
# directly retain high-quality annotations; assign QC-filtered cells by kNN transfer in the reference PCA space
common=ref.var_names
q=raw[:,common].copy()
sc.pp.normalize_total(q,target_sum=1e4); sc.pp.log1p(q)
# reference PCA fitted after per-gene centering/scaling
X=q.X.toarray() if hasattr(q.X,'toarray') else q.X
mu=ref.var.loc[common,'mean'].to_numpy(); sd=ref.var.loc[common,'std'].to_numpy()
X=(X-mu)/sd; np.clip(X,-10,10,out=X)
q_pca=X @ ref.varm['PCs']
clf=KNeighborsClassifier(n_neighbors=15,weights='distance',metric='euclidean')
clf.fit(ref.obsm['X_pca'][:,:40],ref.obs['cell_type'])
all_labels=pd.Series(clf.predict(q_pca[:,:40]),index=raw.obs_names,name='cell_type')
# exact labels for cells retained in QC (avoid any classifier discrepancies)
all_labels.loc[ref.obs_names]=ref.obs['cell_type'].astype(str)
out=pd.DataFrame({'barcode':raw.obs_names.astype(str),'cell_type':all_labels.reindex(raw.obs_names).values})
out.to_csv('labels.csv',index=False)
print('output',out.shape,'missing',out.cell_type.isna().sum(),'barcodes unique',out.barcode.is_unique)
print(out.cell_type.value_counts().to_string())
print('QC retained',ref.n_obs,'transferred',raw.n_obs-ref.n_obs)
print(out.head())