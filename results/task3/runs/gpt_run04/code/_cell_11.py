import pandas as pd, os
x=pd.read_csv('labels.csv')
print(x.head().to_string(index=False))
print('\nFile:',os.path.abspath('labels.csv'),os.path.getsize('labels.csv'),'bytes')
print('columns=',x.columns.tolist(),'n=',len(x),'all nonempty=',x.notna().all().all())
print('\nmarker file:',os.path.getsize('cluster_markers.csv'),'bytes; annotated object:',os.path.getsize('pbmc_annotated.h5ad'),'bytes')