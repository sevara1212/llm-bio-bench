import pandas as pd

df = pd.read_csv('labels.csv')
print("Columns:", df.columns.tolist())
print("Null values:\n", df.isnull().sum())
print("Head:\n", df.head())
print("Tail:\n", df.tail())
print("Shape:", df.shape)