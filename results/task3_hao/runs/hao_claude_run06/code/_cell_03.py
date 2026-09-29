import scanpy as sc
adata = sc.read_h5ad('raw_counts.h5ad')
print(adata.shape)
print(adata.X.max(), adata.X.min())
# check mito genes
mito = [g for g in adata.var_names if g.startswith('MT-')]
print(len(mito), mito[:10])