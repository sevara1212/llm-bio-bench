import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed_adata.h5ad')

# Look at marker gene expression across clusters for res 0.8 and 0.9 or 1.0
marker_genes = ['IL7R', 'CD4', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'MS4A1', 'CD14', 'FCGR3A', 'MS4A7', 'CST3', 'PPBP', 'FCER1A']

print("Marker genes present in raw:")
print([g for g in marker_genes if g in adata.raw.var_names])

# Let's inspect res 0.9 / 1.0 clusters and their marker expression
for res in ['leiden_0.8', 'leiden_0.9', 'leiden_1.0']:
    print(f"\n--- {res} ---")
    sc.tl.rank_genes_groups(adata, res, method='wilcoxon', use_raw=True)
    res_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
    print(res_df)