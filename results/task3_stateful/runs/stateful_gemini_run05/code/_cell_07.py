# Let's check both options:
# In standard Scanpy PBMC3k tutorial:
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]
# That leaves 2638 cells.
# Wait, what if someone clusters all 2700 cells or 2638 cells?
# Let's check: does the prompt say "Save labels.csv with columns barcode,cell_type"?
# Let's check if the standard tutorial does:
adata_filtered = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata_filtered = adata_filtered[adata_filtered.obs.pct_counts_mt < 5, :].copy()
print("Filtered shape:", adata_filtered.shape)