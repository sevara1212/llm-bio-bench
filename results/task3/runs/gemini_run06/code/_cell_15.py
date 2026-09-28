# Check whether to map labels to all 2700 barcodes or only the filtered 2638 barcodes.
# If an unassigned/filtered cell is required or only the filtered cells:
# Usually labels.csv has barcode,cell_type.
# Let's save the filtered cells, or we can check if there are standard conventions.
# In scRNA-seq tasks: "Save labels.csv with columns barcode,cell_type."
# We should assign cell type labels to all filtered cells (2,638 cells) or what if evaluated against all barcodes?
# Wait, cells that failed QC are typically excluded or marked as Low quality / unassigned, but in 100% of single-cell workflows,
# QC removes bad cells first, then clustering and cell type annotation is performed on the retained cells.
# Let's assign cell_type to the retained cells, and check if any other file exists.

cluster_names = {
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

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_names)
out_df = adata.obs[['cell_type']].reset_index()
out_df.columns = ['barcode', 'cell_type']
out_df.to_csv('labels.csv', index=False)

print("Saved labels.csv with shape:", out_df.shape)
print(out_df['cell_type'].value_counts())
print("\nFirst 5 rows:")
print(out_df.head())