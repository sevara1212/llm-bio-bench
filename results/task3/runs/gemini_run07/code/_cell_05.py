# Let's check the prompt instructions carefully:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). 
# Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. 
# Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."

# Notice: "Save labels.csv with columns barcode,cell_type."
# Does every barcode in the input need a label, or only filtered cells?
# Let's check how many barcodes: 2700 barcodes are in filtered_gene_bc_matrices/hg19/barcodes.tsv.
# If we filter cells, some barcodes might be dropped. But wait! If labels.csv is evaluated against ground truth, what does the evaluation expect?
# Does it expect all 2700 barcodes or filtered barcodes?
# Let's check if there are test scripts or other hints.