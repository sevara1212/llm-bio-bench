import scanpy as sc, pandas as pd, numpy as np
from scipy.sparse import diags
ad=sc.read_h5ad('pbmc_qc_clustered.h5ad'); key='leiden_0.5'
labels={'0':'Cycling cytotoxic lymphocyte','1':'CD16 monocyte','2':'Platelet','3':'CD1C dendritic cell','4':'MAIT cell','5':'Erythrocyte','6':'CD4 memory T cell','7':'Memory B cell','8':'Naive CD4 T cell','9':'GZMK+ CD8 T cell','10':'GZMB+ plasmacytoid dendritic cell','11':'Naive B cell','12':'CLEC9A dendritic cell','13':'Erythroid progenitor','14':'Classical monocyte','15':'Cytotoxic T cell','16':'Activated gamma-delta T cell','17':'Gamma-delta T cell','18':'Plasma cell','19':'Natural killer cell','20':'Plasmacytoid dendritic cell'}
ad.obs['cell_type']=ad.obs[key].map(labels).astype(str)
# properly normalized raw state and marker ranking
counts=ad.layers['counts']; tot=np.asarray(counts.sum(1)).ravel(); x=diags(1e4/tot).dot(counts); x.data=np.log1p(x.data); ad.X=x; ad.raw=ad.copy()
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=100,pts=True)
markers=[]
for cl in ad.obs[key].cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=cl,key='rank_genes_groups').head(20);d.insert(0,'cell_type',labels[cl]);d.insert(0,'cluster',cl);markers.append(d)
pd.concat(markers).to_csv('cluster_markers.csv',index=False)
raw=sc.read_h5ad('raw_counts.h5ad');out=pd.DataFrame({'barcode':raw.obs_names.astype(str)});out['cell_type']=out.barcode.map(ad.obs['cell_type']).astype(object)
assign={'E2L8_CCCGGAATCGGCTGTG':'Platelet','E2L4_GAGTGAGCAAACCATC':'Low-quality cell','E2L1_CTGAATGAGCGGGTTA':'Low-quality cell','E2L2_GGAACCCGTATGTCTG':'Low-quality cell'}
out['cell_type']=out['cell_type'].fillna(out['barcode'].map(assign)); assert out.cell_type.notna().all() and len(out)==7841
out.to_csv('labels.csv',index=False);ad.write_h5ad('pbmc_qc_clustered_annotated.h5ad')
print('PASS',out.shape, 'unique barcodes',out.barcode.nunique(), 'null',out.cell_type.isna().sum())
print(out.cell_type.value_counts().to_string())
print('\nmarkers exists',pd.read_csv('cluster_markers.csv').shape)