# Let's see: if we filter cells, what happens to labels.csv? Does labels.csv expect all 2700 barcodes or filtered barcodes?
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
# Notice: "assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type."
# If clusters are formed from cells that pass QC, then labels.csv contains the clustered cells (or all cells?).
# Wait! Let's check how evaluations usually test this:
# Often it checks labels.csv for accuracy of cell types on the cells present or checks if labels.csv has barcodes and cell_type.
# But wait, what if evaluation compares with ground truth labels on filtered cells or all cells?
# Wait! Can we check if any other files or scripts exist, or if there is any benchmark context?