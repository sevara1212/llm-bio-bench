# Let's inspect those high hb cells
print("high hb cells (>20%):", (adata.obs['pct_counts_hb'] > 20).sum())
print("high hb cells (>5%):", (adata.obs['pct_counts_hb'] > 5).sum())