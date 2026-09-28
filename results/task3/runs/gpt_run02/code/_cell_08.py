import scanpy as sc, pandas as pd, numpy as np
# Recover all barcodes in original input; annotate retained Leiden clusters
raw=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols'); raw.var_names_make_unique()
proc=sc.read_h5ad('pbmc_processed.h5ad')
labels={'0':'CD4 T cell','1':'CD14+ monocyte','2':'B cell','3':'CD8 T cell','4':'FCGR3A+ monocyte','5':'NK cell','6':'dendritic cell','7':'platelet'}
# all source barcodes, explicit QC label for six excluded barcodes
out=pd.DataFrame({'barcode':raw.obs_names.astype(str)})
maplab=proc.obs['leiden_0.8'].astype(str).map(labels)
out['cell_type']=out.barcode.map(maplab).fillna('Low-quality cell')
out.to_csv('labels.csv',index=False)
# cluster audit
summary=(proc.obs['leiden_0.8'].astype(str).value_counts().rename_axis('cluster').reset_index(name='n_cells'))
summary['cell_type']=summary.cluster.map(labels)
summary=summary.sort_values('cluster',key=lambda s:s.astype(int))
summary.to_csv('cluster_annotations.csv',index=False)
print(summary.to_string(index=False))
print('\nlabels shape:',out.shape,'unique barcodes:',out.barcode.nunique(),'missing labels:',out.cell_type.isna().sum())
print(out.cell_type.value_counts().to_string())
print('\nCSV preview:');print(pd.read_csv('labels.csv').head().to_string(index=False))