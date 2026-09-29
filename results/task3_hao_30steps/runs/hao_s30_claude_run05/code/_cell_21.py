import scanpy as sc
adata = sc.read_h5ad('full_with_clusters.h5ad')
genes_check = ['FOXP3','IL2RA','CTLA4','GATA3']
genes_check = [g for g in genes_check if g in adata.var_names]
import scanpy as sc
df = sc.get.obs_df(adata, keys=genes_check+['leiden'])
print(df.groupby('leiden').mean())