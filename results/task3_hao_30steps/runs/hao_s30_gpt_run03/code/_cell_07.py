import scanpy as sc
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.6';sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=30)
r=ad.uns['rank_genes_groups']
with open('markers_by_cluster.txt','w') as f:
 for cl in ad.obs[key].cat.categories:
  s=f'{cl}\t(n={sum(ad.obs[key]==cl)})\t'+', '.join(r['names'][cl][:25])+'\n';f.write(s);print(s,end='')