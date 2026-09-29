# Let's inspect all clusters in detail!
for col in markers.columns:
    print(f"Cluster {col}: {list(markers[col])}")