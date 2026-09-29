import pandas as pd
import scanpy as sc

raw = sc.read_h5ad('raw_counts.h5ad')
labels = pd.read_csv('labels.csv')

print("Raw cells:", raw.n_obs)
print("Labels rows:", labels.shape[0])
print("Unique barcodes in labels:", labels.barcode.nunique())
print("Any duplicated barcodes:", labels.barcode.duplicated().sum())
print("All labels barcodes exist in raw:", labels.barcode.isin(raw.obs_names).all())
print(labels['cell_type'].value_counts())