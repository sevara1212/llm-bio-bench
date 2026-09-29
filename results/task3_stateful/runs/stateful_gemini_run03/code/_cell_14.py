# Let's check markers for cluster 0, 1, 4, 5
# Cluster 0 has CD3D, LDHB -> CD4 T cells (or T cells)
# Let's check CD4, CD8A, CD8B in cluster 0 vs 1
marker_check = ['CD3D', 'CD4', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'MS4A1', 'CD14', 'FCGR3A', 'CST3', 'PPBP', 'IL7R']
for gene in marker_check:
    if gene in adata.raw.var_names:
        print(f"Gene {gene} in raw")

# Let's check expression by cluster:
import scanpy as sc
df_expr = pd.DataFrame(index=adata.obs['leiden'].unique())
for gene in marker_check:
    if gene in adata.raw.var_names:
        adata.obs[gene] = adata.raw[:, gene].X.toarray().flatten()
        print(adata.obs.groupby('leiden')[gene].mean())