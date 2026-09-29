# Let's run a standard pipeline on adata without filtering out cells so every barcode is preserved,
# OR check if standard PBMC preprocessing does filtering.
# Notice: "Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
# If we filter cells, what happens to their labels in labels.csv?
# If labels.csv has all 7841 barcodes, it's safest to either keep all 7841 or mark low quality.
# But let's check standard PBMC clustering first.

# Save raw counts in a layer
adata.layers['counts'] = adata.X.copy()
# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # keep full normalized/log1p in raw for marker genes

# Feature selection
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", adata.var['highly_variable'].sum())

# PCA
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_0.5')
sc.tl.leiden(adata, resolution=0.8, key_added='leiden_0.8')

print("Clusters at res 0.5:", adata.obs['leiden_0.5'].value_counts())
print("Clusters at res 0.8:", adata.obs['leiden_0.8'].value_counts())