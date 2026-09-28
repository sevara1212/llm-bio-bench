import pandas as pd, os
x=pd.read_csv('labels.csv')
print(x.head(10).to_string(index=False))
print('unique barcodes',x.barcode.nunique(),'rows',len(x),'duplicate barcodes',x.barcode.duplicated().sum())
print('file bytes',os.path.getsize('labels.csv'))
print(pd.crosstab(x.cell_type,columns='cells').to_string())