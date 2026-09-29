# Let's check cluster 0 further: can it be split into Naive CD4 and Memory CD4 T cells?
# In Seurat pbmc3k tutorial:
# 0: Naive CD4 T
# 1: CD14+ Mono
# 2: Memory CD4 T
# 3: B
# 4: CD8 T
# 5: FCGR3A+ Mono
# 6: NK
# 7: DC
# 8: Platelet
# Let's check resolution=1.2 or 1.4 to see if cluster 0 splits
sc.tl.leiden(adata, resolution=1.3)
print(adata.obs['leiden'].value_counts())