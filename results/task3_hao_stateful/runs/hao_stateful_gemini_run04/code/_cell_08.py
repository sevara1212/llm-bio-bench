# Let's see if this dataset comes from a public PBMC benchmark or tutorial
# Check some gene names and barcodes
print("Gene names sample:", adata.var_names[:20])
# Check batches
adata.obs['batch'] = [x.split('_')[0] for x in adata.obs.index]
print(adata.obs['batch'].value_counts())