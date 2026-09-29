# Let's inspect barcodes and potential metadata
print(adata.obs.index[:10])
# Look for potential existing scripts, check what _cell_01.py or _kernel.py is
with open('_cell_01.py') as f:
    print("cell_01:", f.read()[:500])