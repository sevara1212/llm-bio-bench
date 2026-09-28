import scanpy as sc
import pandas as pd
import numpy as np

# Load original raw data
raw_adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
raw_adata.var_names_make_unique()
raw_barcodes = raw_adata.obs_names.tolist()

# Load processed adata
adata = sc.read_h5ad('adata_processed.h5ad')
sc.tl.leiden(adata, resolution=1.0, key_added='leiden')

# Check CCR7 and S100A4
for g in ['CCR7', 'S100A4']:
    e0 = adata.raw[:, g].X.toarray().flatten()[adata.obs['leiden'] == '0'].mean()
    e6 = adata.raw[:, g].X.toarray().flatten()[adata.obs['leiden'] == '6'].mean()
    print(f"{g}: cluster 0 = {e0:.3f}, cluster 6 = {e6:.3f}")

# Cell type mapping
# 0: CD4 T cells
# 1: CD14+ Monocytes
# 2: B cells
# 3: CD8 T cells
# 4: FCGR3A+ Monocytes
# 5: NK cells
# 6: CD4 T cells (or CD4+ T cells)
# 7: Dendritic cells
# 8: Megakaryocytes

mapping = {
    '0': 'CD4 T cells',
    '1': 'CD14+ Monocytes',
    '2': 'B cells',
    '3': 'CD8 T cells',
    '4': 'FCGR3A+ Monocytes',
    '5': 'NK cells',
    '6': 'CD4 T cells',
    '7': 'Dendritic cells',
    '8': 'Megakaryocytes'
}

adata.obs['cell_type'] = adata.obs['leiden'].map(mapping)

# Build the labels dataframe
# Note: barcodes filtered out during QC
filtered_out_barcodes = [bc for bc in raw_barcodes if bc not in adata.obs_names]
print(f"Retained cells: {len(adata.obs_names)}, Filtered out cells: {len(filtered_out_barcodes)}")

# Let's map all barcodes
# For filtered out cells, we can either label them or keep only filtered, but let's check standard practice.
# Usually labels.csv contains the filtered cells or all cells. If all cells, filtered out could be "Low quality" or "Filtered".
# But if it says "Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type."
# We should save barcodes and their cell type labels.
df_labels = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type']
})

df_labels.to_csv('labels.csv', index=False)
print("labels.csv created with shape:", df_labels.shape)
print(df_labels.head())
print(df_labels['cell_type'].value_counts())