# Let's check the prompt instructions carefully:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). 
# Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. 
# Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
#
# Note: "Save labels.csv with columns barcode,cell_type."
# Does every barcode in the input need a label, or only the filtered ones?
# Let's check how many cells are in raw_counts: 7841.
# Notice min n_genes is 499, max pct_counts_mt is 15.06 (only 3 cells > 15%).
# Could this dataset ALREADY have been lightly pre-filtered, or is standard filtering expected?
# If we filter cells out, their barcodes wouldn't be in labels.csv, OR should all 7841 be in labels.csv?
# Let's check standard PBMC QC:
# e.g., min_genes=200 (all cells have >= 499, so 0 removed if min_genes=200!)
# min_cells=3 for genes
# What if pct_counts_mt < 15 or 20? If < 20, 0 cells removed! If min_genes=200, 0 cells removed!
print(adata.n_obs)