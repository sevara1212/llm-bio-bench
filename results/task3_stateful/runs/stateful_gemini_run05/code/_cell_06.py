# Let's inspect genes and cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print("After min_genes=200 and min_cells=3:", adata.shape)

# Let's see if we filter out mt > 5% or keep standard scanpy tutorial filtering:
# In scanpy tutorial:
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]
# Wait, if we filter out cells, do we need to assign a cell type label to every cluster or every cell?
# "assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type."
# If a barcode was filtered out, what should its label be? Or does the evaluation expect only the retained cells or all 2700 cells?