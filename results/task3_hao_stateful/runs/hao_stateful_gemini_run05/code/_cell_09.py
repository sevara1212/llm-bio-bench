for col in markers_03.columns:
    print(f"Cluster {col}: {', '.join(markers_03[col].values[:7])}")