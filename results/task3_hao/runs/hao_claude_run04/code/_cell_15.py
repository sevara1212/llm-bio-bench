import pandas as pd
labels = pd.read_csv('labels.csv')
print(labels.shape)
print(labels['barcode'].duplicated().sum(), "duplicate barcodes")
print(labels['cell_type'].isna().sum(), "missing cell types")
print(labels['cell_type'].value_counts())