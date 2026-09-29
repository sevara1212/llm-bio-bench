# Let's run a pipeline on adata without filtering cells first, to see what clusters form!
# If we filter cells, will labels.csv need to have ALL 7841 barcodes or only filtered ones?
# Notice: "Save labels.csv with columns barcode,cell_type."
# If we filter cells, those cells wouldn't have a cluster/cell_type!
# Let's check the wording:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
# Let's see if 7841 cells are already pre-filtered quality cells!
# min genes is 499 (basically 500!), min counts is 898 (~900), max mt is 15%.
# That means 7841 barcodes are ALREADY the filtered cell barcodes from CellRanger or standard filtering!