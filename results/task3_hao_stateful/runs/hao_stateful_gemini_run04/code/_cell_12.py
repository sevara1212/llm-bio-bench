# Let's inspect all clusters from 0 to 18
for col in marker_df_05.columns:
    print(f"Cluster {col}: {', '.join(marker_df_05[col].tolist())}")