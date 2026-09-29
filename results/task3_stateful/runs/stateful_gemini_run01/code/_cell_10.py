# Let's check marker genes in the tutorial and in our clusters at res=1.0 vs res=0.8 vs res=0.5
# In the classic scanpy pbmc3k tutorial:
# Marker genes:
# IL7R, CCR7 -> CD4+ T
# CD14, LYZ -> CD14+ Monocytes
# MS4A1 -> B cells
# CD8A -> CD8+ T
# FCGR3A, MS4A7 -> FCGR3A+ Monocytes
# GNLY, NKG7 -> NK cells
# FCER1A, CST3 -> Dendritic cells
# PPBP -> Megakaryocytes

# Let's inspect res=0.8 and res=1.0 with rank_genes_groups
adata_processed.raw = adata
sc.tl.leiden(adata_processed, resolution=1.0, random_state=0)
sc.tl.rank_genes_groups(adata_processed, 'leiden', method='wilcoxon')
result = adata_processed.uns['rank_genes_groups']
groups = result['names'].dtype.names
marker_df = pd.DataFrame(
    {group + '_' + key: result[key][group][:5]
    for group in groups for key in ['names']}
)
print("Resolution 1.0 markers:")
print(marker_df)