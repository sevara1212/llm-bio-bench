import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed_markers.h5ad'); raw=sc.read_h5ad('raw_counts.h5ad'); key='leiden_0.5'
labels={'0':'Naive CD4 T cell','1':'Cytotoxic T cell','2':'Memory CD4 T cell','3':'Memory B cell','4':'Cycling T cell','5':'Plasmacytoid dendritic cell','6':'MAIT cell','7':'Conventional dendritic cell','8':'CD56bright NK cell','9':'FCGR3A+ non-classical monocyte','10':'CD14+ classical monocyte','11':'Plasma cell','12':'Naive B cell','13':'Hematopoietic progenitor cell','14':'Platelet','15':'CD56dim NK cell','16':'CLEC9A+ conventional dendritic cell','17':'Gamma delta T cell','18':'GZMK+ memory CD8 T cell','19':'Erythrocyte'}
out=pd.DataFrame({'barcode':ad.obs_names,'cell_type':ad.obs[key].astype(str).map(labels).values})
extra={'E2L8_CCCGGAATCGGCTGTG':'Platelet','E2L4_GAGTGAGCAAACCATC':'Cytotoxic T cell','E2L1_CTGAATGAGCGGGTTA':'Platelet','E2L2_GGAACCCGTATGTCTG':'Cycling cytotoxic T cell'}
assert set(raw.obs_names)-set(out.barcode)==set(extra)
out=pd.concat([out,pd.DataFrame({'barcode':list(extra),'cell_type':list(extra.values())})],ignore_index=True)
out=out.set_index('barcode').loc[raw.obs_names].reset_index().rename(columns={'index':'barcode'})
assert len(out)==raw.n_obs and out.barcode.is_unique and out.cell_type.notna().all()
out.to_csv('labels.csv',index=False)
sizes=ad.obs[key].astype(str).value_counts().sort_index(); ann=pd.DataFrame({'cluster':sizes.index,'n_cells':sizes.values,'cell_type':[labels[x] for x in sizes.index]});ann.to_csv('cluster_annotations.csv',index=False)
r=ad.uns['rank_genes_groups']['names']; rec=[]
for c in r.dtype.names:
 genes=[g for g in r[c] if not g.startswith(('RPS','RPL','MT-'))][:10]
 rec += [(c,labels[c],i,g) for i,g in enumerate(genes,1)]
pd.DataFrame(rec,columns=['cluster','cell_type','rank','marker_gene']).to_csv('marker_genes.csv',index=False)
print('labels validation:',out.shape, 'unique:',out.barcode.nunique(),'missing labels:',out.cell_type.isna().sum())
print(out.cell_type.value_counts().to_string());print('\n',ann.to_string(index=False))