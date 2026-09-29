# Let's check Ding et al. Nature Biotech 2020 pbmc data or similar benchmarking dataset?
# Ding et al. systematic comparison of single-cell and single-nucleus RNA-seq methods:
# pbmc1, pbmc2? E2L1-8 (e.g. 10x v2 / v3 / etc. or 10x Chromium lanes?).
# Let's check gene filtering:
# Standard QC filtering:
# Often cells are filtered, but wait! The instruction says:
# "Save labels.csv with columns barcode,cell_type."
# Does labels.csv expect all 7841 barcodes or only filtered ones?
# Usually, evaluation scripts check labels.csv against all barcodes in the dataset (or index on adata.obs_names)!
# If we filter cells out, those cells won't have a label in labels.csv, or if the benchmark compares on the input adata, missing barcodes could cause key errors or lower accuracy!