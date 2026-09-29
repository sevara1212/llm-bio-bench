import scanpy as sc

adata = sc.read_h5ad('processed.h5ad')
adata.raw = None
# Set adata.raw to normalized logged counts without ribo/mito if we want or just check top logfoldchanges
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon', use_raw=False)

markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'LEF1', 'IL7R', 'S100A4', 
           'FOXP3', 'IL2RA', 'TRDC', 'TRGC1', 'NKG7', 'GNLY', 'NCAM1', 'FCGR3A',
           'MS4A1', 'CD19', 'CD79A', 'IGHM', 'IGHD', 'CD27', 'MZB1', 'SDC1',
           'CD14', 'FCGR3A', 'LYZ', 'S100A9', 'MS4A7', 'FCER1A', 'CST3', 'CLEC9A', 'CLEC4C', 'IL3RA', 'LILRA4',
           'PPBP', 'PF4', 'HBB', 'HBA1', 'MKI67', 'TOP2A', 'CD34']

# Print expression of these markers across clusters
import pandas as pd
df_expr = pd.DataFrame(index=range(21))
for m in markers:
    if m in adata.var_names:
        df_expr[m] = [adata[adata.obs['leiden_0.5'] == str(i), m].X.mean() for i in range(21)]

print(df_expr[['CD3D', 'CD4', 'CD8A', 'CCR7', 'IL7R', 'CCL5' if 'CCL5' in adata.var_names else 'CD3E']].round(2))