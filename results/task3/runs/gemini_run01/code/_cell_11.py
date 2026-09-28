import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed_adata.h5ad')

# Compare cluster 0 and cluster 1 in leiden_0.9:
# Which one is CD4 T cells and which is what?
# Let's check marker genes between cluster 0 and 1, or markers of cluster 0 vs rest, cluster 1 vs rest:
sc.tl.rank_genes_groups(adata, 'leiden_0.9', groups=['0', '1'], reference='rest', method='wilcoxon', use_raw=True)
print("Top 10 for cluster 0:")
print(adata.uns['rank_genes_groups']['names']['0'][:10])
print("\nTop 10 for cluster 1:")
print(adata.uns['rank_genes_groups']['names']['1'][:10])

# Also check specific markers: CCR7, S100A4, SELL, FOXP3, IL7R, CD8A, CD8B
t_markers = ['IL7R', 'CCR7', 'S100A4', 'SELL', 'CD4', 'CD8A', 'CD8B', 'CCL5', 'GZMK', 'GZMB', 'GZMA', 'NKG7']
raw_df = pd.DataFrame(adata.raw[:, [g for g in t_markers if g in adata.raw.var_names]].X.toarray(), 
                      columns=[g for g in t_markers if g in adata.raw.var_names], 
                      index=adata.obs_names)
raw_df['cluster'] = adata.obs['leiden_0.9'].values
print("\nT cell markers in clusters 0, 1, 4, 6:")
print(raw_df[raw_df['cluster'].isin(['0', '1', '4', '6'])].groupby('cluster').mean())