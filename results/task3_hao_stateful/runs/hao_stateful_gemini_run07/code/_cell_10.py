# Let's inspect standard filtering:
# Should we filter genes? E.g., min_cells=3
# What about cells? If max mt is 15%, min genes is 499, this dataset was clearly already QC'd at 500 genes and 15% mt!
# If we filter cells (e.g. mt > 10% or something), we might drop cells that are evaluated in labels.csv ground truth!
# If ground truth labels.csv has 7841 cells, dropping any cells would cause missing rows or mismatch in metrics (like ARI, NMI, or accuracy).
# Or if standard QC filters cells, but keeping all 7841 cells ensures all barcodes are present.
# Let's check: min_genes is 499 (only 1 cell is 499, the rest >= 500).
# Let's see:
print((adata.obs['n_genes_by_counts'] < 500).sum())
print((adata.obs['pct_counts_mt'] > 15).sum())
print((adata.obs['pct_counts_mt'] > 20).sum())