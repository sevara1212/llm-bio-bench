import scanpy as sc, pandas as pd, os
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.7'
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=30)
r=ad.uns['rank_genes_groups']
for c in ['6','22']:
 gs=[str(x) for x in r['names'][c] if not str(x).startswith(('RPL','RPS','MT-'))][:20]
 print(c,(ad.obs[key]==c).sum(),', '.join(gs))
# File verification
x=pd.read_csv('labels.csv'); raw=sc.read_h5ad('raw_counts.h5ad')
print('labels file',os.path.getsize('labels.csv'),'bytes; schema',x.columns.tolist(),'unique rows',x.barcode.nunique())
print('barcode sets match',set(x.barcode)==set(raw.obs_names),'null labels',x.cell_type.isna().sum())
print(pd.read_csv('cluster_markers.csv').groupby(['cluster','cell_type']).size().to_string())