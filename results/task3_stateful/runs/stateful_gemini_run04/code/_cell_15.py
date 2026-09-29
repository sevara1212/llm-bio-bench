# Save to labels.csv
df_labels.to_csv('labels.csv', index=False)
print("labels.csv saved successfully!")
print("Checking file exists and first lines:")
print(pd.read_csv('labels.csv').head(10))