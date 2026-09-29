import scanpy as sc

ad = sc.read_h5ad('clustered.h5ad')

sc.tl.leiden(ad, resolution=0.5, key_added='leiden_05')
print(ad.obs['leiden_05'].value_counts())

ad.write('clustered.h5ad')