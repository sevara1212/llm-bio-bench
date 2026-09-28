# Let's inspect the canonical PBMC tutorial in Scanpy!
# In scanpy tutorial:
# marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
#                 'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
#                 'FCGR3A', 'MS4A7', 'CST3', 'FCER1A', 'PPBP']
# Usually in Scanpy pbmc3k tutorial:
# leiden resolution 1.0 or similar has 8 clusters:
# 0: CD4 T cells
# 1: CD14+ Monocytes
# 2: B cells
# 3: CD8 T cells
# 4: NK cells
# 5: FCGR3A+ Monocytes
# 6: Dendritic Cells
# 7: Megakaryocytes

# Let's check higher resolutions:
for res in [0.7, 0.8, 0.9, 1.0, 1.2]:
    sc.tl.leiden(adata_sub, resolution=res, key_added=f'leiden_{res}', random_state=0)
    print(f"Res {res}: {adata_sub.obs[f'leiden_{res}'].nunique()} clusters")