import pandas as pd
import scanpy as sc

# Final checks
labels = pd.read_csv('labels.csv')
raw = sc.read_h5ad('raw_counts.h5ad')

print("Total barcodes in labels.csv:", len(labels))
print("Total barcodes in raw data:", raw.n_obs)
print("Match:", set(labels['barcode']) == set(raw.obs_names))
print("No duplicates:", labels['barcode'].duplicated().sum() == 0)
print("\nUnique cell types:")
print(labels['cell_type'].value_counts())
print("\nAny missing values:", labels.isna().sum().sum())