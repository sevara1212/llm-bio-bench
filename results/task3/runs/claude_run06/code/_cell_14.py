import pandas as pd

labels = pd.read_csv('labels.csv')
print("Shape:", labels.shape)
print("Columns:", labels.columns.tolist())
print("Any duplicated barcodes?", labels['barcode'].duplicated().sum())
print("Any missing values?", labels.isnull().sum().sum())
print()
print(labels['cell_type'].value_counts())
print()
print("Sample rows:")
print(labels.sample(5, random_state=42))