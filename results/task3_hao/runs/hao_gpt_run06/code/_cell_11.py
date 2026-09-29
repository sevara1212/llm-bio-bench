import scanpy as sc, numpy as np, pandas as pd
from scipy.sparse import diags
# load clustered (with retained cells) and restore correctly normalized expression as .raw for marker testing
ad=sc.read_h5ad('pbmc_qc_clustered.h5ad')
counts=ad.layers['counts'].copy(); totals=np.asarray(counts.sum(axis=1)).ravel()
lognorm=diags(1e4/totals).dot(counts); lognorm.data=np.log1p(lognorm.data)
ad.X=lognorm
ad.raw=ad.copy()
# cluster labels, selected resolution 0.5 after inspecting markers
key='leiden_0.5'
labels={
'0':'Cycling cytotoxic lymphocyte',
'1':'CD16 monocyte',
'2':'Platelet',
'3':'CD1C dendritic cell',
'4':'MAIT cell',
'5':'Erythrocyte',
'6':'CD4 memory T cell',
'7':'Memory B cell',
'8':'Naive CD4 T cell',
'9':'GZMK+ CD8 T cell',
'10':'GZMB+ plasmacytoid dendritic cell',
'11':'Naive B cell',
'12':'CLEC9A dendritic cell',
'13':'Erythroid progenitor',
'14':'Classical monocyte',
'15':'Cytotoxic T cell',
'16':'Activated gamma-delta T cell',
'17':'Gamma-delta T cell',
'18':'Plasma cell',
'19':'Natural killer cell',
'20':'Plasmacytoid dendritic cell'}
ad.obs['cell_type']=ad.obs[key].map(labels).astype('category')
# DE uses correctly normalized raw
sc.tl.rank_genes_groups(ad,groupby=key,method='wilcoxon',use_raw=True,n_genes=100,pts=True)
rows=[]
for cl in ad.obs[key].cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=cl,key='rank_genes_groups')
 d.insert(0,'cell_type',labels[cl]);d.insert(0,'cluster',cl)
 rows.append(d.head(20))
markers=pd.concat(rows,ignore_index=True)
markers.to_csv('cluster_markers.csv',index=False)
# labels for retained + failed QC raw cells. Give three damaged cells explicit low-quality label; platelet exception.
raw=sc.read_h5ad('raw_counts.h5ad')
out=pd.DataFrame({'barcode':raw.obs_names.astype(str)})
out['cell_type']=out.barcode.map(ad.obs['cell_type'])
# determine remaining in raw and assign based on top genes / QC review
missing=out.cell_type.isna()
print('missing',out.loc[missing,'barcode'].tolist())
assign={'E2L8_CCCGGAATCGGCTGTG':'Platelet',
        'E2L4_GAGTGAGCAAACCATC':'Low-quality cell',
        'E2L1_CTGAATGAGCGGGTTA':'Low-quality cell',
        'E2L2_GGAACCCGTATGTCTG':'Low-quality cell'}
out.loc[missing,'cell_type']=out.loc[missing,'barcode'].map(assign)
assert out.cell_type.notna().all() and len(out)==raw.n_obs and out.barcode.is_unique
out.to_csv('labels.csv',index=False)
# summary and persist final object
ad.write_h5ad('pbmc_qc_clustered_annotated.h5ad')
print(out.cell_type.value_counts().to_string())
print('labels rows',len(out),'markers rows',len(markers))
print(markers[['cluster','cell_type','names','scores','pvals_adj']].groupby(['cluster','cell_type']).head(5).to_string(index=False))