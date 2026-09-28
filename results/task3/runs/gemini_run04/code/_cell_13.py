# Check CD14, FCGR3A, MS4A7 in cluster 2 vs 4
for gene in ['CD14', 'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'CD3D', 'CD4', 'CD8A', 'GNLY', 'NKG7', 'MS4A1']:
    print(gene, "present?", gene in adata.raw.var_names)

# Let's check mean expression in clusters for res 0.8
import pandas as pd
genes_to_check = ['CD14', 'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'CD3D', 'CD4', 'CD8A', 'GNLY', 'NKG7', 'MS4A1', 'PPBP']
expr = pd.DataFrame(adata.raw[:, genes_to_check].X.toarray(), columns=genes_to_check, index=adata.obs_names)
expr['cluster'] = adata.obs['leiden_0.8']
print("\nMean expression by cluster (res 0.8):")
print(expr.groupby('cluster').mean().round(3))