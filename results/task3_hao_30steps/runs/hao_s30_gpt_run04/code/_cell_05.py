import scanpy as sc, pandas as pd
adata=sc.read_h5ad('qc_processed.h5ad')
sc.tl.rank_genes_groups(adata,groupby='leiden_0.5',method='wilcoxon',use_raw=True,n_genes=12)
r=adata.uns['rank_genes_groups']
for c in adata.obs['leiden_0.5'].cat.categories:
 print('\n',c,adata.obs['leiden_0.5'].value_counts()[c], list(r['names'][c][:12]), list(r['scores'][c][:5].round(1)))