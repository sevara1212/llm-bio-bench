import scanpy as sc
adata = sc.read_h5ad('processed.h5ad')
print("Leiden 0.5 clusters counts:")
print(adata.obs['leiden_0.5'].value_counts())
print("\nLeiden 0.8 clusters counts:")
print(adata.obs['leiden_0.8'].value_counts())