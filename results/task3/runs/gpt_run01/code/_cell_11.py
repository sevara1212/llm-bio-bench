import scanpy as sc, numpy as np, pandas as pd
# Load original and derive QC; use retained analysis object for labels
raw=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols'); raw.var_names_make_unique()
raw.var['mt']=raw.var_names.str.startswith('MT-'); sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True,log1p=False)
ad=sc.read_h5ad('pbmc_processed.h5ad')
# resolution 1.0 labels; manual annotations based on marker genes
mapping={
 '0':'CD4 T cell',       # IL7R, LTB, CD3D
 '1':'Naive CD4 T cell', # CCR7, LTB, T-cell genes
 '2':'Classical monocyte', # LYZ, S100A8/9, FCN1
 '3':'CD8 T cell',       # CD3D, CCL5, GZMK/NKG7
 '4':'B cell',           # MS4A1, CD79A/B, CD74
 '5':'Non-classical monocyte', # LST1, FCGR3A, LILRB/CTSS
 '6':'Natural killer cell', # NKG7, GNLY, PRF1, GZMB
 '7':'Proliferating lymphocyte', # STMN1, HMGB2, T-cell genes
 '8':'Platelet'          # PPBP, PF4, GNG11
}
cl=ad.obs['leiden_1.0'].astype(str)
out=pd.DataFrame({'barcode':ad.obs_names, 'cell_type':cl.map(mapping)})
# add the two QC-excluded cells by 15-NN classification in PCA space using shared normalized/log data and stored PC loading
excluded=raw.obs_names.difference(ad.obs_names)
# calculate PCs from analysis loadings after matching genes; transform normalized log matrix
rr=raw[:,ad.var_names].copy();sc.pp.normalize_total(rr,target_sum=1e4);sc.pp.log1p(rr)
pcs=np.asarray(rr.X @ ad.varm['PCs'][:,:40])
from sklearn.neighbors import KNeighborsClassifier
train=np.array([raw.obs_names.get_loc(x) for x in ad.obs_names])
test=np.array([raw.obs_names.get_loc(x) for x in excluded])
knn=KNeighborsClassifier(n_neighbors=15,weights='distance').fit(pcs[train],out.cell_type)
pred=knn.predict(pcs[test])
extra=pd.DataFrame({'barcode':excluded,'cell_type':pred})
out=pd.concat([out,extra],ignore_index=True)
# preserve original barcode order; checks
out['barcode']=pd.Categorical(out.barcode,categories=raw.obs_names,ordered=True);out=out.sort_values('barcode');out['barcode']=out['barcode'].astype(str)
out.to_csv('labels.csv',index=False)
print('Output:',out.shape,' unique barcode:',out.barcode.nunique(),' missing labels:',out.cell_type.isna().sum())
print(out.cell_type.value_counts().to_string())
print('QC-excluded/NN assigned:');print(extra.to_string(index=False))
print('file checks:',pd.read_csv('labels.csv').head(),pd.read_csv('labels.csv').columns.tolist())