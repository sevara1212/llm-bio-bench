# Keep track of original barcodes
all_barcodes = pd.Index(adata.obs_names)
print("Total original barcodes:", len(all_barcodes))