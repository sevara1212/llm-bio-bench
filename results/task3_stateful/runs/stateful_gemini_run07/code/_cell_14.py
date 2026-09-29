# Let's see: resolution 0.5 gave 6 clusters. What if resolution 0.6 or 0.8?
# Let's check with resolution 0.8 (which is standard scanpy PBMC tutorial!):
# In standard scanpy tutorial:
# sc.tl.leiden(adata) # default resolution is 1.0 or 0.8? In older scanpy default is 1.0. In tutorial:
# Actually, the tutorial uses resolution: default is 1.0!
# Wait, let's see what clusters appear at resolution 1.0 or 0.5.
# At resolution 0.5:
# C0: CD3D+, CD8A/B low/moderate, CD4 -> T cells (CD4+ T cells / memory T cells)
# C1: CD14+, FCGR3A+, CST3+ -> Monocytes
# C2: NKG7+, GNLY+, CD8A+, CD8B+, CD3D+ -> NK cells / Cytotoxic T cells (NK and CD8 T cells)
# C3: MS4A1+, CD79A+ -> B cells
# C4: CST3+, FCER1A+, HLA-DRA+ -> Dendritic cells
# C5: PPBP+, PF4+ -> Megakaryocytes / Platelets

# Wait! What if we check resolution 1.0?
sc.tl.leiden(adata, resolution=0.8, random_state=0)
print("Res 0.8 clusters:", adata.obs['leiden'].nunique())
sc.tl.leiden(adata, resolution=1.0, random_state=0)
print("Res 1.0 clusters:", adata.obs['leiden'].nunique())