import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed_pbmc.h5ad'); key='leiden_0.6'
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',n_genes=30)
r=ad.uns['rank_genes_groups']
rows=[]
for cl in r['names'].dtype.names:
 rows.append([cl,ad.obs[key].value_counts()[cl],', '.join(r['names'][cl][:20])])
pd.DataFrame(rows,columns=['cluster','n','top_markers']).to_csv('cluster_markers.csv',index=False)
print(pd.DataFrame(rows,columns=['cluster','n','top_markers']).to_string(index=False))