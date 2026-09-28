import pandas as pd
labels = pd.read_csv('labels.csv')
print(labels.shape)
print(labels['cell_type'].unique())
print(labels.isnull().sum())
print(labels['barcode'].duplicated().sum())