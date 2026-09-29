import scanpy as sc
adata = sc.read_h5ad('full_with_clusters.h5ad')

# Check cluster 9 vs 23 -- both look pDC-like from IL3RA/LILRA4/CLEC4C, let's compare more carefully
genes_check = ['IL3RA','LILRA4','CLEC4C','GZMB','MZB1','JCHAIN','TCF4','FCER1A','CD79A','STMN1']
df = sc.get.obs_df(adata, keys=genes_check+['leiden'])
sub = df[df['leiden'].isin(['9','23'])]
print(sub.groupby('leiden').mean())
print()
print(adata.obs['leiden'].value_counts()[['9','23']])