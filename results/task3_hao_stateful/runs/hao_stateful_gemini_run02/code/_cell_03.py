# Note the prompt: "Save labels.csv with columns barcode,cell_type"
# Wait, are we supposed to filter cells or assign a label to EVERY barcode in raw_counts.h5ad or filtered cells?
# Let's re-read carefully:
# "This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood mononuclear cells (PBMCs). Perform standard quality control, normalisation and clustering, identify marker genes, and assign a cell type label to every cluster. Save labels.csv with columns barcode,cell_type. Work step by step and check your outputs."
# Wait, if we filter cells, what happens to filtered cells? Or is it already pre-filtered to 7841 cells?
# Let's check the QC metrics of the 7841 cells.
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
print(adata.obs.describe())