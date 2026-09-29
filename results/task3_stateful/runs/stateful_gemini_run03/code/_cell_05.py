# Check standard scanpy tutorial filtering:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

# But wait, what if labels.csv needs to contain ALL barcodes or only filtered barcodes?
# Let's re-read the prompt:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
# Notice: the directory is "filtered_gene_bc_matrices", which already filtered empty droplets (2700 cells is the classic PBMC 3k dataset!).
# In PBMC 3k: 2700 cells.
# If cells are filtered out during QC, what should happen to their barcodes in labels.csv? Or should QC filter outliers, or keep all cells?
# Let's check how many cells are filtered with standard thresholds:
print("Cells with n_genes < 200:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("Cells with n_genes >= 2500:", (adata.obs['n_genes_by_counts'] >= 2500).sum())
print("Cells with pct_counts_mt >= 5:", (adata.obs['pct_counts_mt'] >= 5).sum())