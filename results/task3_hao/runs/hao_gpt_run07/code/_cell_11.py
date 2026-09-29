import pandas as pd, os
x=pd.read_csv('labels.csv')
print(x.head().to_string(index=False)); print(x.columns.tolist(),len(x),x.barcode.duplicated().sum(),x.cell_type.isna().sum())
print('files:',[(f,os.path.getsize(f)) for f in ['labels.csv','cluster_markers.csv','cluster_annotation_report.csv','processed_pbmc.h5ad']])