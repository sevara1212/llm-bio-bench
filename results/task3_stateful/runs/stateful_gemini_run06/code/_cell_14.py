df_out = pd.DataFrame({
    'barcode': adata_scaled.obs_names,
    'cell_type': adata_scaled.obs['cell_type']
})

df_out.to_csv('labels.csv', index=False)
print("Saved labels.csv with shape:", df_out.shape)
print("Preview:")
print(pd.read_csv('labels.csv').head(10))