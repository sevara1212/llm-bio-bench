# Let's inspect the instructions carefully:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). 
# Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. 
# Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."

# Notice: "assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type."
# If we filter cells during QC, what should labels.csv contain?
# Does labels.csv require all barcodes from raw_counts.h5ad, or only filtered cells, or should QC filtering be conservative or assign filtered cells?
# Let's check how many cells are in raw_counts.h5ad: 7841 cells.
# Let's see the distributions of n_genes, total_counts, pct_counts_mt.