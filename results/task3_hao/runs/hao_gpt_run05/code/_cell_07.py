import scanpy as sc
ad=sc.read_h5ad('processed.h5ad');r=ad.uns['rank_genes_groups']; key='leiden_0.4'
for c in ad.obs[key].cat.categories:
 print(c, sum(ad.obs[key]==c), '|', ','.join(r['names'][c][:8]))