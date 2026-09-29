import scanpy as sc
import anndata as ad

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("obs head:", adata.obs.head())
print("var head:", adata.var.head())
if os.path.exists('_cell_01.py'):
    with open('_cell_01.py') as f:
        print("cell_01:", f.read())