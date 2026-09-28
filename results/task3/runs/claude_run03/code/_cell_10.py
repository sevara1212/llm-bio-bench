import pandas as pd

labels = pd.read_csv('labels.csv')
print(labels.shape)
print(labels['cell_type'].value_counts())
print(labels.head(10))
print("Unique barcodes:", labels['barcode'].nunique())
print("Any NaN cell_type:", labels['cell_type'].isna().sum())