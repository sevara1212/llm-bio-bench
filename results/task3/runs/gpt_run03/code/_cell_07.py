import scanpy as sc, pandas as pd, numpy as np
from sklearn.neighbors import KNeighborsClassifier
from scipy import sparse
# processed clusters and annotation
A=sc.read_h5ad('pbmc_processed_r09.h5ad')
label_map={'0':'CD4 T cell','1':'CD14 monocyte','2':'B cell','3':'Naive CD4 T cell','4':'CD8 T cell','5':'FCGR3A+ monocyte','6':'NK cell','7':'Dendritic cell'}
A.obs['cell_type']=A.obs['leiden'].map(label_map)
# Save ranked marker results, including statistics
markers=[]
for c in sorted(A.obs.leiden.unique(),key=int):
 d=sc.get.rank_genes_groups_df(A,group=c).head(20)
 d.insert(0,'cell_type',label_map[c]); d.insert(0,'cluster',c); markers.append(d)
pd.concat(markers,ignore_index=True).to_csv('cluster_markers.csv',index=False)
# Project 26 QC-excluded cells to the annotated population using kNN on normalized log1p HVG expression.
raw=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols'); raw.var_names_make_unique()
raw.var['mt']=raw.var_names.str.upper().str.startswith('MT-'); sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True)
keep=(raw.obs.n_genes_by_counts>=300)&(raw.obs.n_genes_by_counts<=2500)&(raw.obs.pct_counts_mt<10)
# Verify barcode matching/order consistency
assert set(raw.obs_names[keep])==set(A.obs_names)
# same library-size normalization/log transform applied before selecting reference HVGs
sc.pp.normalize_total(raw,target_sum=1e4);sc.pp.log1p(raw)
hv=A.var_names[A.var.highly_variable]
Xref=raw[A.obs_names,hv].X
Xall=raw[:,hv].X
knn=KNeighborsClassifier(n_neighbors=15,weights='distance',metric='cosine',n_jobs=1)
knn.fit(Xref,A.obs.cell_type)
pred=knn.predict(Xall)
# Retained labels use direct cluster annotation, never projected labels
lab=pd.Series(pred,index=raw.obs_names);lab.loc[A.obs_names]=A.obs.cell_type
out=pd.DataFrame({'barcode':raw.obs_names,'cell_type':lab.loc[raw.obs_names].values})
out.to_csv('labels.csv',index=False)
print('Output:',out.shape, 'unique barcodes',out.barcode.nunique(),'missing labels',out.cell_type.isna().sum())
print(out.cell_type.value_counts().to_string())
print('QC-excluded=',(~keep).sum(),' projected labels:',lab.loc[~keep].value_counts().to_dict())
print('\nmarkers file rows:',len(pd.concat(markers)),' labels header/sample:')
print(pd.read_csv('labels.csv').head().to_string(index=False))
print('\nCluster annotation:')
print(pd.DataFrame({'cluster':list(label_map),'cell_type':list(label_map.values()),'n_cells':[int((A.obs.leiden==x).sum()) for x in label_map]}).to_string(index=False))