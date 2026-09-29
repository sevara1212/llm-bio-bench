import pickle,pandas as pd,scanpy as sc
with open('clustered.pkl','rb') as f: ad=pickle.load(f)
labels={
'0':'CD4 naive T cell','1':'CD8 effector memory T cell','2':'CD4 memory T cell','3':'Memory B cell','4':'Cycling lymphocyte','5':'CD4 memory T cell','6':'CD16 monocyte','7':'Conventional dendritic cell 2','8':'Natural killer cell','9':'Plasmacytoid dendritic cell','10':'Plasma cell','11':'Naive B cell','12':'Hematopoietic progenitor cell','13':'CD14 monocyte','14':'Platelet','15':'Natural killer cell','16':'Conventional dendritic cell 1','17':'Gamma delta T cell','18':'CD8 memory T cell','19':'Erythroid cell','20':'Plasmacytoid dendritic cell'}
ad.obs['cell_type']=ad.obs.cluster.map(labels).astype(str)
# output QC-pass labels
out=pd.DataFrame({'barcode':ad.obs_names,'cell_type':ad.obs.cell_type})
# include filtered-out barcodes as QC-failed to ensure each input barcode has a result
raw=sc.read_h5ad('raw_counts.h5ad'); raw.var['mt']=raw.var_names.str.upper().str.startswith('MT-');sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True)
full=pd.DataFrame({'barcode':raw.obs_names})
full=full.merge(out,on='barcode',how='left');full['cell_type']=full.cell_type.fillna('Low-quality cell (mitochondrial RNA >15%)')
full.to_csv('labels.csv',index=False)
# All DE results concise table
markers=[]
for cl in sorted(labels,key=int):
 d=sc.get.rank_genes_groups_df(ad,group=cl).head(20).copy();d.insert(0,'cell_type',labels[cl]);d.insert(0,'cluster',cl);markers.append(d)
pd.concat(markers).to_csv('cluster_markers.csv',index=False)
print(full.shape,full.cell_type.value_counts().to_string())
print('unique',full.barcode.nunique(),'missing',full.cell_type.isna().sum())
print('filtered',full[full.cell_type.str.startswith('Low')].to_string(index=False))
print('markers',pd.read_csv('cluster_markers.csv').shape)