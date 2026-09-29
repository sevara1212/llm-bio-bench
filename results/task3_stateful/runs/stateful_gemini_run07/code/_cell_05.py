# Check how many cells are filtered with standard cutoffs
# Standard PBMC3k tutorial:
# n_genes < 2500 (or pct_counts_mt < 5% in standard tutorial, adata = adata[adata.obs.pct_counts_mt < 5, :])
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# But wait! Does the prompt require labels.csv for ALL barcodes or only filtered barcodes, or how is it usually structured?
# Prompt: "Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
# Let's check if there are any instructions or existing scripts