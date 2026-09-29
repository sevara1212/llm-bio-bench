# Let's map the clusters to cell types:
# 0: CD4 T cells
# 1: CD14+ Monocytes
# 2: B cells
# 3: CD8 T cells
# 4: FCGR3A+ Monocytes
# 5: NK cells
# 6: CD4 T cells (or T cells)
# 7: Dendritic cells
# 8: Megakaryocytes

# Let's check cluster_name mapping:
cluster_map = {
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

adata_filtered.obs['cell_type'] = adata_filtered.obs['leiden'].map(cluster_map)

# Now what about cells that were filtered out if labels.csv needs all 2700 barcodes or just filtered cells?
# Let's save both or check how labels.csv is usually expected.
# "Save labels.csv with columns barcode,cell_type."
# Usually in single-cell benchmarking / pipelines, labels.csv contains the analyzed cells:
df_labels = pd.DataFrame({
    'barcode': adata_filtered.obs_names,
    'cell_type': adata_filtered.obs['cell_type'].values
})
df_labels.to_csv('labels.csv', index=False)
print("Saved labels.csv with shape:", df_labels.shape)
print(df_labels.head(10))
print(df_labels['cell_type'].value_counts())