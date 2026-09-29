# Let's check cluster 6 vs cluster 0
# Both are CD3D+, IL7R+.
# In standard PBMC 3k:
# 0: CD4 T cells
# 1: CD14+ Monocytes
# 2: B cells
# 3: CD8 T cells
# 4: FCGR3A+ Monocytes
# 5: NK cells
# 6: Dendritic cells or CD4 T cells? Wait!
# Let's check what cluster 6 expresses:
sc.tl.rank_genes_groups(adata_filtered, 'leiden', groups=['6'], reference='0', method='wilcoxon')
print("6 vs 0:", adata_filtered.uns['rank_genes_groups']['names']['6'][:10])
# What about naive vs memory CD4 T?
# In scanpy tutorial:
# new_cluster_names = [
#     'CD4 T', 'CD14+ Monocytes',
#     'B', 'CD8 T',
#     'FCGR3A+ Monocytes', 'NK',
#     'Dendritic', 'Megakaryocytes'
# ]
# Wait, in the tutorial there were 8 clusters! Here we have 9 clusters (0 to 8).
# Cluster 7 has FCER1A (2.15), HLA-DRA, CST3 -> Dendritic Cells!
# Cluster 8 has PPBP (5.69) -> Megakaryocytes!
# Cluster 6 has high ribosomal genes! Let's check CD4 vs CD8 or naive CD4 T cells.
# Let's check CCR7, SELL (CD62L), S100A4 in cluster 0 vs 6:
for gene in ['CCR7', 'SELL', 'S100A4', 'IL7R', 'FOXP3']:
    if gene in adata_filtered.raw.var_names:
        idx = adata_filtered.raw.var_names.get_loc(gene)
        val0 = np.mean(adata_filtered.raw.X[(adata_filtered.obs['leiden'] == '0').values, idx].toarray())
        val6 = np.mean(adata_filtered.raw.X[(adata_filtered.obs['leiden'] == '6').values, idx].toarray())
        print(f"{gene}: cl0={val0:.2f}, cl6={val6:.2f}")