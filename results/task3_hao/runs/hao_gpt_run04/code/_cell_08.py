import pickle,scanpy as sc,pandas as pd,numpy as np
from sklearn.metrics.pairwise import cosine_similarity
with open('analysis_results.pkl','rb') as f:a=pickle.load(f)
k='leiden_0.6'
labels={'0':'Naive CD4 T cell','1':'Memory CD4 T cell','2':'Cytotoxic CD8 T cell','3':'Naive B cell','4':'CD161+ memory T cell','5':'CD1c+ conventional dendritic cell','6':'Memory B cell','7':'FCGR3A+ non-classical monocyte','8':'Natural killer cell','9':'Plasmacytoid dendritic cell','10':'Erythroid progenitor','11':'CD14+ classical monocyte','12':'Plasma cell','13':'Natural killer cell','14':'Cycling natural killer cell','15':'Platelet','16':'Cycling lymphocyte','17':'CLEC9A+ conventional dendritic cell','18':'Gamma delta T cell','19':'GZMK+ memory CD8 T cell','20':'Erythrocyte'}
a.obs['cell_type']=a.obs[k].map(labels).astype(str)
print(a.obs.cell_type.value_counts().to_string())
# Map the 146 filtered cells to retained clusters using cosine similarity of log-normalized expression over HVGs.
raw=sc.read_h5ad('raw_counts.h5ad')
raw.var['mt']=raw.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True,log1p=False)
keep=(raw.obs.n_genes_by_counts>=500)&(raw.obs.pct_counts_mt<12)
assert raw.obs_names[keep].equals(a.obs_names)
removed=raw[~keep].copy(); sc.pp.normalize_total(removed,target_sum=1e4);sc.pp.log1p(removed)
# only HVGs used for robust mapping
idx=np.where(a.var.highly_variable.values)[0]
cent=[]; cats=a.obs[k].cat.categories
for c in cats: cent.append(np.asarray(a[a.obs[k]==c].X[:,idx].mean(axis=0)).ravel())
sims=cosine_similarity(removed.X[:,idx],np.vstack(cent))
assigned=np.array(cats)[sims.argmax(1)]
print('removed',removed.n_obs,'mapping',pd.Series(assigned).value_counts().to_dict(),'mean max cos',sims.max(1).mean().round(3),'min',sims.max(1).min().round(3))
# every original barcode, ordered original input
out=pd.DataFrame({'barcode':raw.obs_names.astype(str),'cell_type':pd.Series(index=raw.obs_names,dtype='object')})
out.loc[a.obs_names,'cell_type']=a.obs['cell_type'].values
out.loc[removed.obs_names,'cell_type']=[labels[x] for x in assigned]
assert out.cell_type.notna().all() and out.barcode.is_unique and len(out)==raw.n_obs
out.to_csv('labels.csv',index=False)
print(out.head().to_string(index=False));print('output',out.shape, 'types',out.cell_type.nunique())