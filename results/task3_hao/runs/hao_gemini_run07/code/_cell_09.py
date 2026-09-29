import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Let's see if we filter lightly or keep all
# If an automated test checks labels.csv against all barcodes, len(labels) == 7841 is safest or filtered?
# Wait! Let's check: "Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
# Does the task say: "Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type."

# Let's inspect what clusters we get without dropping cells vs dropping cells.
adata_proc = adata.copy()
# Save raw counts in a layer or raw
adata_proc.raw = adata_proc

sc.pp.normalize_total(adata_proc, target_sum=1e4)
sc.pp.log1p(adata_proc)
sc.pp.highly_variable_genes(adata_proc, n_top_genes=2000, subset=False)
adata_hvg = adata_proc[:, adata_proc.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata_hvg, random_state=42)
sc.tl.leiden(adata_hvg, resolution=0.5, random_state=42)

print(adata_hvg.obs['leiden'].value_counts())