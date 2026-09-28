# Let's check how both 0 and 1 compare to Naive CD4 vs Memory CD4:
# In Seurat PBMC3k tutorial:
# 0: Naive CD4+ T (IL7R, CCR7)
# 1: CD14+ Mono (LYZ, CD14)
# 2: Memory CD4+ (IL7R, S100A4)
# 3: B (MS4A1)
# 4: CD8+ T (CD8A)
# 5: FCGR3A+ Mono (FCGR3A, MS4A7)
# 6: NK (GNLY, NKG7)
# 7: DC (FCER1A, CST3)
# 8: Platelet (PPBP)

# In our leiden_0.9:
# Cluster 0: IL7R+, S100A4 high, LTB high -> Memory CD4+ T cell
# Cluster 1: IL7R+, CCR7 high, S100A4 lower -> Naive CD4+ T cell
# Cluster 2: CD14+, LYZ+, S100A9+ -> CD14+ Monocyte
# Cluster 3: MS4A1+, CD79A+, CD79B+ -> B cell
# Cluster 4: CD8A+, CD8B+, CCL5+ -> CD8+ T cell
# Cluster 5: FCGR3A+, MS4A7+, LST1+, CST3+ -> FCGR3A+ Monocyte (Non-classical monocyte)
# Cluster 6: GNLY+, NKG7+, PRF1+, GZMB+ -> NK cell
# Cluster 7: HLA-DRA+, HLA-DPA1+, FCER1A+, CST3+ -> Dendritic cell
# Cluster 8: PPBP+, PF4+ -> Megakaryocyte / Platelet

# Let's verify cluster sizes in leiden_0.9:
import scanpy as sc
adata = sc.read_h5ad('processed_adata.h5ad')
print(adata.obs['leiden_0.9'].value_counts())