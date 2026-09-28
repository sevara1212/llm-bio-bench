import pandas as pd
labels = pd.read_csv('labels.csv')
print(labels.shape)
print(labels.isnull().sum())
print(labels['barcode'].duplicated().sum())
print(labels['cell_type'].unique())
print(labels.head(10))