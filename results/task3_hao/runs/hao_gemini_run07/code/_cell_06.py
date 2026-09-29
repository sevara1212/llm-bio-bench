import scanpy as sc
import pandas as pd

raw = sc.read_h5ad('raw_counts.h5ad')
print("Total cells in raw:", raw.n_obs)
print("Any batch or sample in index?", raw.obs_names[:10])
# Look at whether the task expects all 7841 barcodes in labels.csv:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type."