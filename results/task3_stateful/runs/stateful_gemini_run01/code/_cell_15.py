# Let's map cluster labels to cell types:
# 0: CD4+ T cells
# 1: CD14+ Monocytes
# 2: B cells
# 3: CD8+ T cells
# 4: FCGR3A+ Monocytes
# 5: NK cells
# 6: CD4+ T cells
# 7: Dendritic cells
# 8: Megakaryocytes

cluster_to_cell_type = {
    '0': 'CD4+ T cells',
    '1': 'CD14+ Monocytes',
    '2': 'B cells',
    '3': 'CD8+ T cells',
    '4': 'FCGR3A+ Monocytes',
    '5': 'NK cells',
    '6': 'CD4+ T cells',
    '7': 'Dendritic cells',
    '8': 'Megakaryocytes'
}

adata_processed.obs['cell_type'] = adata_processed.obs['leiden'].map(cluster_to_cell_type)
print(adata_processed.obs['cell_type'].value_counts())

# Now let's check barcode coverage:
# adata_processed has 2638 cells.
# What about the 62 filtered cells (from original 2700)?
# Usually labels.csv contains the analyzed cells or all barcodes. Let's inspect raw barcodes:
raw_barcodes = pd.read_csv('filtered_gene_bc_matrices/hg19/barcodes.tsv', header=None, sep='\t')[0].values
print("Total raw barcodes:", len(raw_barcodes))

# Check how labels.csv is formatted
labels_df = pd.DataFrame({
    'barcode': adata_processed.obs_names,
    'cell_type': adata_processed.obs['cell_type'].values
})
labels_df.to_csv('labels.csv', index=False)
print("Saved labels.csv with shape:", labels_df.shape)
print(labels_df.head(10))