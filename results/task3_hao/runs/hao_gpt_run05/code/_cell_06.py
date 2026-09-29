import scanpy as sc
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.4'
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=15)
r=ad.uns['rank_genes_groups']
for c in ad.obs[key].cat.categories:
 print('\nCLUSTER',c,'n=',sum(ad.obs[key]==c))
 print(', '.join(r['names'][c][:15]))
ad.write_h5ad('processed.h5ad')