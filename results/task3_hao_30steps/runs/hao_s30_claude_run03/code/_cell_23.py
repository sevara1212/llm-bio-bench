import scanpy as sc
import pandas as pd

ad = sc.read_h5ad('final_labeled.h5ad')
raw = sc.read_h5ad('raw_counts.h5ad')

# Build labels df for cells passing QC
labels_df = pd.DataFrame({
    'barcode': ad.obs_names,
    'cell_type': ad.obs['cell_type'].values
})

# Identify filtered-out barcodes and label them
all_barcodes = set(raw.obs_names)
passed_barcodes = set(ad.obs_names)
filtered_out = all_barcodes - passed_barcodes
print("Filtered out barcodes:", filtered_out)

filtered_df = pd.DataFrame({
    'barcode': list(filtered_out),
    'cell_type': ['Low quality (filtered out)'] * len(filtered_out)
})

final_labels = pd.concat([labels_df, filtered_df], ignore_index=True)
print(final_labels.shape)
print(final_labels['cell_type'].value_counts())

# Reorder to match raw barcode order for cleanliness
final_labels = final_labels.set_index('barcode').loc[raw.obs_names].reset_index()
print(final_labels.head())
print(final_labels.shape)

final_labels.to_csv('labels.csv', index=False)
print("Saved labels.csv")