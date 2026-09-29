# Note: The prompt asks:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). 
# Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. 
# Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."

# Wait, if we filter out cells during QC, what happens to labels.csv? Does labels.csv need to contain EVERY barcode in raw_counts.h5ad, or only filtered cells?
# "assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type."
# Let's check how many barcodes there are in raw_counts.h5ad: 7841.
# Notice min n_genes is 499, max mt is 15.06%. The data may ALREADY be pre-filtered for low quality cells!
# Look at min n_genes = 499, min total counts = 898, max mt = 15.06%.
# This strongly suggests the data has already been pre-filtered for low-quality cells / empty droplets (e.g. n_genes >= 500, mt < 15%)!