import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed.h5ad')
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')

canonical_markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'MS4A1', 'CD79A', 'GNLY', 'NKG7', 'NCAM1',
                     'CD14', 'FCGR3A', 'CST3', 'CLEC4C', 'LILRA4', 'PPBP', 'PF4', 'FCER1A', 'CD1C']
available_markers = [m for m in canonical_markers if m in adata.raw.var_names]
print("Canonical markers available:", available_markers)

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

for g in groups:
    genes = result['names'][g][:8]
    print(f"Cluster {g} (n={sum(adata.obs['leiden_0.5'] == g)}): {', '.join(genes)}")