import pandas as pd

labels = pd.read_csv('labels.csv')
print(labels.shape)
print(labels['cell_type'].value_counts())
print("Any duplicated barcodes:", labels['barcode'].duplicated().sum())
print("Any missing values:", labels.isnull().sum())
print(labels.head())