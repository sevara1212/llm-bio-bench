import scanpy as sc, pandas as pd, numpy as np
# Load clustering and associate annotations
ad=sc.read_h5ad('pbmc_processed.h5ad')
cluster_to_type={
 '0':'Memory CD4 T cell',
 '1':'CD14+ monocyte',
 '2':'Naive CD4 T cell',
 '3':'B cell',
 '4':'Cytotoxic T cell',
 '5':'FCGR3A+ monocyte',
 '6':'Natural killer cell',
 '7':'Dendritic cell',
 '8':'Platelet',
 '9':'Cycling T/NK cell'
}
ad.obs['cell_type']=ad.obs['leiden_1.0'].astype(str).map(cluster_to_type)
print(ad.obs.groupby(['leiden_1.0','cell_type'],observed=True).size())
# Create labels for every source barcode, explicitly retain QC failures as QC label
raw=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols',cache=False)
labels=pd.DataFrame({'barcode':raw.obs_names, 'cell_type':'Low-quality cell'})
labels.loc[labels.barcode.isin(ad.obs_names),'cell_type']=labels.loc[labels.barcode.isin(ad.obs_names),'barcode'].map(ad.obs['cell_type'])
labels.to_csv('labels.csv',index=False)
print('\nlabels:',labels.shape, 'missing',labels.cell_type.isna().sum(),'duplicates',labels.barcode.duplicated().sum())
print(labels.cell_type.value_counts())
# Save marker table from log-normalized, non-scaled expression (Wilcoxon)
mk=ad.copy()
# rebuild log normalized X directly from QC-filtered raw (using same barcode/gene ordering)
rawsub=raw[ad.obs_names,ad.var_names].copy()
sc.pp.normalize_total(rawsub,target_sum=1e4); sc.pp.log1p(rawsub)
rawsub.obs['cluster']=ad.obs['leiden_1.0'].values
sc.tl.rank_genes_groups(rawsub,groupby='cluster',method='wilcoxon',pts=True,key_added='markers')
markers=sc.get.rank_genes_groups_df(rawsub,group=None,key='markers')
markers['cell_type']=markers['group'].astype(str).map(cluster_to_type)
markers.to_csv('cluster_markers.csv',index=False)
print('\nTop markers (first 8)')
for cl in sorted(cluster_to_type,key=int):
 print(cl,cluster_to_type[cl],':',', '.join(markers[markers.group==cl].head(8)['names']))
# Persist analysis annotation
ad.write('pbmc_annotated.h5ad')