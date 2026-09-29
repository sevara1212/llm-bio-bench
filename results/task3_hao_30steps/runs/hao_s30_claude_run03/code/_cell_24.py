import pandas as pd

final_labels = pd.read_csv('labels.csv')
print(final_labels.columns.tolist())
final_labels = final_labels.rename(columns={'index': 'barcode'})
final_labels.to_csv('labels.csv', index=False)
print(final_labels.head())
print(final_labels.shape)
print(final_labels.isna().sum())