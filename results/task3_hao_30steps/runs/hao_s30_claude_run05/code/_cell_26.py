import pandas as pd

labels = pd.read_csv('labels.csv')
print("Shape:", labels.shape)
print("Columns:", labels.columns.tolist())
print("Unique barcodes:", labels['barcode'].nunique())
print("Any NaN:", labels.isna().sum().sum())
print()
print("Cell type distribution:")
print(labels['cell_type'].value_counts())
print()
print("Sample rows:")
print(labels.head(10))

# Verify barcode set matches original data
import scanpy as sc
adata_orig = sc.read_h5ad('raw_counts.h5ad')
print()
print("Original n_obs:", adata_orig.n_obs)
print("Labels n_obs:", labels.shape[0])
print("Difference (filtered out cells):", adata_orig.n_obs - labels.shape[0])