# Let's inspect the prompt carefully:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs).
# Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster.
# Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
# Notice: "Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
# Should labels.csv contain ALL 7841 barcodes, or only filtered barcodes, or what is standard?
# Wait! Let's check if standard QC filters cells, or if filtering cells means barcodes not in filtered set are excluded OR if QC is just checking/filtering.
# Wait, let's see why total_counts max is 67792, min genes 499, max mt is 15.06%!
# Notice min n_genes is 499, max mt is 15.06%. This dataset might ALREADY have been lightly pre-filtered (e.g. min genes 500 or 200, mt < 15%)!
# Let's check:
print("Min genes:", adata.obs['n_genes_by_counts'].min())
print("Max mt:", adata.obs['pct_counts_mt'].max())