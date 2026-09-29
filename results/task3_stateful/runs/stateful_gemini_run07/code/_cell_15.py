# Let's inspect res 1.0 (9 clusters is EXACTLY the Scanpy tutorial clusters: 8 or 9 clusters!)
# In Scanpy PBMC3k tutorial:
# new_cluster_names = [
#     'CD4 T', 'CD14 Monocytes', 'B', 'CD8 T', 
#     'NK', 'FCGR3A Monocytes', 'Dendritic', 'Megakaryocytes'
# ]
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markers_df = pd.DataFrame({group: [result['names'][group][i] for i in range(8)] for group in sorted(groups, key=int)})
print(markers_df)

# Let's check markers across all 9 clusters:
key_genes = ['IL7R', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP', 'CD3D']
for g in key_genes:
    if g in adata.raw.var_names:
        idx = adata.raw.var_names.get_loc(g)
        vals = adata.raw.X[:, idx].toarray().flatten()
        print(f"{g:8s}:", [f"C{c}:{vals[adata.obs['leiden'].values == str(c)].mean():.1f}" for c in range(9)])

# Let's assign cell types based on these:
# Let's see:
# Cluster 0: IL7R high, CD3D high -> CD4+ T cells
# Cluster 1: CD14 high, LYZ high -> CD14+ Monocytes
# Cluster 2: CD79A, MS4A1 high -> B cells
# Cluster 3: CD8A high, CD3D high -> CD8+ T cells
# Cluster 4: GNLY, NKG7 high, CD3D low/absent -> NK cells
# Cluster 5: FCGR3A, MS4A7, CST3 high -> FCGR3A+ Monocytes
# Cluster 6: FCER1A, CST3, HLA high -> Dendritic cells
# Cluster 7: PPBP, PF4 high -> Megakaryocytes
# Cluster 8: Let's see what cluster 8 is!