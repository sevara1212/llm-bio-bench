# Let's check the origin of this dataset: barcodes start with L3_, E2L7_, E2L2_, E2L4_, L5_
# This looks like a multiplexed or multi-sample PBMC dataset (e.g. 10x Genomics PBMC or Kang et al. or Hao et al. Azimuth PBMC or similar)!
# Let's check some gene names:
print("CD3D in var:", 'CD3D' in adata.var_names)
print("MS4A1 in var:", 'MS4A1' in adata.var_names)
print("CD14 in var:", 'CD14' in adata.var_names)
print("FCGR3A in var:", 'FCGR3A' in adata.var_names)
print("GNLY in var:", 'GNLY' in adata.var_names)
print("PPBP in var:", 'PPBP' in adata.var_names)