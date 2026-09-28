import scanpy as sc, pandas as pd, numpy as np, os
# Re-read raw data to make barcode-level output and exactly reproduce QC flags
raw=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols'); raw.var_names_make_unique()
raw.var['mt']=raw.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True,percent_top=None,log1p=False)
pass_qc=(raw.obs.n_genes_by_counts>=200)&(raw.obs.n_genes_by_counts<=2500)&(raw.obs.pct_counts_mt<5)
# Cluster identities generated after QC + library-size normalisation/log1p/HVG/PCA/neighbors/Leiden (resolution 0.8)
ad=sc.read_h5ad('pbmc_processed_hvg.h5ad')
cluster=ad.obs['leiden_0.8'].astype(str)
celltypes={'0':'CD4 T cell','1':'CD14 monocyte','2':'B cell','3':'CD8 cytotoxic T cell','4':'FCGR3A+ monocyte','5':'Natural killer cell','6':'Naive CD4 T cell','7':'Dendritic cell','8':'Platelet'}
labels=pd.Series('Low-quality cell',index=raw.obs_names,dtype='object')
labels.loc[ad.obs_names]=cluster.map(celltypes).values
out=pd.DataFrame({'barcode':raw.obs_names,'cell_type':labels.values})
out.to_csv('labels.csv',index=False)
# Marker table (from a fresh rank test)
ad.obs['cluster']=cluster
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=True,n_genes=50)
markers=sc.get.rank_genes_groups_df(ad,group=None)
markers['cell_type']=markers['group'].map(celltypes)
markers.to_csv('cluster_markers_wilcoxon.csv',index=False)
print('labels exists',os.path.exists('labels.csv'),'shape',out.shape,'unique barcode',out.barcode.nunique())
print(out.cell_type.value_counts().to_string())
print('QC passed',pass_qc.sum(),'failed',(~pass_qc).sum(),'labels not low-quality',(out.cell_type!='Low-quality cell').sum())
print('nulls',out.isna().sum().to_dict())
print('\nSelected top nontechnical markers:')
for c in map(str,range(9)):
 m=markers[markers.group==c].names.tolist()
 m=[x for x in m if not x.startswith(('RPL','RPS','MT-')) and x not in ['MALAT1','TPT1','LDHB','EEF1A1','B2M']]
 print(c,celltypes[c],':',', '.join(m[:8]))