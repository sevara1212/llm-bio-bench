# Map clusters to cell types
cluster_annotations = {
    '0': 'T cells',
    '1': 'NK cells',
    '2': 'CD14+ Monocytes',
    '3': 'B cells',
    '4': 'FCGR3A+ Monocytes',
    '5': 'Dendritic cells',
    '6': 'Megakaryocytes'
}

adata_filtered.obs['cell_type'] = adata_filtered.obs['leiden'].map(cluster_annotations)

# Now check barcodes:
# If user expects all barcodes or filtered barcodes:
# "Save labels.csv with columns barcode,cell_type."
# Usually this is for the quality-filtered cells (2638 cells) or all cells?
# Let's check:
df_labels = pd.DataFrame({
    'barcode': adata_filtered.obs_names,
    'cell_type': adata_filtered.obs['cell_type']
})

print("Filtered shape:", df_labels.shape)
print(df_labels.head())
print(df_labels['cell_type'].value_counts())