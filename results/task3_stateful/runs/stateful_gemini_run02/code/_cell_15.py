# Let's run res=1.0 clustering and map the cluster names
sc.tl.leiden(adata, resolution=1.0)
adata.obs['cell_type'] = adata.obs['leiden'].map({
    '0': 'CD4 T cells',
    '1': 'CD14+ Monocytes',
    '2': 'B cells',
    '3': 'CD8+ T cells',
    '4': 'FCGR3A+ Monocytes',
    '5': 'NK cells',
    '6': 'Dendritic cells',
    '7': 'Platelets'
})

df_out = pd.DataFrame({
    'barcode': adata.obs_names,
    'cell_type': adata.obs['cell_type']
})

df_out.to_csv('labels.csv', index=False)
print("Saved labels.csv with shape:", df_out.shape)
print(df_out['cell_type'].value_counts())
print(df_out.head())