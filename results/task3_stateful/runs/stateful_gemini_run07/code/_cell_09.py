# Let's re-read the exact prompt:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."

# In standard single-cell analysis, QC filters out low-quality cells (or barcodes).
# When you cluster and assign cell type to every cluster, the resulting cells are the filtered cells.
# Wait, or does it assign labels to the clustered cells?
# "assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type."
# Since clusters are defined on the filtered cells, only clustered cells have a cluster and thus a cluster-level cell type label!
# Wait, what if someone evaluated on either:
# Could filtered-out cells be labeled "Filtered" or omitted?
# In Seurat and Scanpy tutorials, filtered cells are removed from the dataset. Then clustering is performed, and clusters are annotated.
# Then `adata.obs['cell_type']` has index (barcode) and value (cell_type).
# Exporting `adata.obs[['cell_type']].to_csv('labels.csv')` or with index_label='barcode'.