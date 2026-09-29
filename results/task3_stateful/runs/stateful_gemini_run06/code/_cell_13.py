cluster_labels = {
    '0': 'Memory CD4 T',
    '1': 'Naive CD4 T',
    '2': 'CD14+ Monocytes',
    '3': 'B cells',
    '4': 'CD8 T',
    '5': 'FCGR3A+ Monocytes',
    '6': 'NK cells',
    '7': 'Dendritic cells',
    '8': 'Megakaryocytes'
}

adata_scaled.obs['cell_type'] = adata_scaled.obs['leiden'].map(cluster_labels)
print(adata_scaled.obs['cell_type'].value_counts())

# Now check barcode output
# Does the task require saving labels.csv for filtered cells or all 2700 cells?
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs).
# Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster.
# Save labels.csv with columns barcode,cell_type."
print("Number of cells in adata_scaled:", len(adata_scaled))