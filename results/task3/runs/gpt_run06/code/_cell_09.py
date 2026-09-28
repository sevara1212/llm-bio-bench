import scanpy as sc, pandas as pd, numpy as np
# derive DE markers from unmodified UMI counts, aligning retained observations to clustering
p='filtered_gene_bc_matrices/hg19/'
raw=sc.read_10x_mtx(p,var_names='gene_symbols',cache=False); raw.var_names_make_unique()
raw.var['mt']=raw.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True)
keep=(raw.obs.n_genes_by_counts>=200)&(raw.obs.pct_counts_mt<10)
cl=sc.read_h5ad('pbmc_qc_clustered.h5ad')
raw=raw[cl.obs_names].copy() # exact retained order
raw.obs['leiden']=cl.obs['leiden'].astype(str).values
sc.pp.normalize_total(raw,target_sum=1e4); sc.pp.log1p(raw)
sc.tl.rank_genes_groups(raw,'leiden',method='wilcoxon',n_genes=25)
rows=[]
for c in sorted(raw.obs.leiden.unique(),key=int):
 genes=list(raw.uns['rank_genes_groups']['names'][c][:15])
 rows.append({'cluster':c,'n_cells':int((raw.obs.leiden==c).sum()),'markers':'; '.join(genes)})
pd.DataFrame(rows).to_csv('cluster_markers.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
# labels for retained cluster annotations; all barcodes, QC-failed labeled Low quality
mapping={'0':'CD4 T cell','1':'NK cell','2':'CD14 monocyte','3':'B cell','4':'FCGR3A+ monocyte','5':'Dendritic cell','6':'Megakaryocyte'}
out=pd.DataFrame({'barcode':cl.obs_names.astype(str),'cell_type':cl.obs.leiden.astype(str).map(mapping).values})
failed=raw.obs_names # irrelevant
allbc=sc.read_10x_mtx(p,var_names='gene_symbols',cache=False).obs_names.astype(str)
out=out.set_index('barcode').reindex(allbc)
out['cell_type']=out.cell_type.fillna('Low quality')
out.reset_index().rename(columns={'index':'barcode'}).to_csv('labels.csv',index=False)
print('\nLabels:',out.cell_type.value_counts().to_dict())
print('output shape',pd.read_csv('labels.csv').shape,'missing',pd.read_csv('labels.csv').isna().sum().sum())