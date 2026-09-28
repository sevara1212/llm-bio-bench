import scanpy as sc, pandas as pd, numpy as np, os
# Load preliminary object: contains 2,698 QC-passing cells and Leiden resolution 1.0
A=sc.read_h5ad('pbmc_processed_prelim.h5ad')
cluster_to_type={
 '0':'CD4 T cell',       # CD3D/E, IL7R, LTB, CCR7
 '1':'CD4 T cell',       # CD3D/E, IL7R, LTB
 '2':'CD14 monocyte',    # LYZ, S100A8/A9, FCN1
 '3':'B cell',           # MS4A1, CD79A/B, CD74
 '4':'CD8 T cell',       # CD3D, CD8A, CCL5, NKG7
 '5':'FCGR3A+ monocyte', # LST1, FCGR3A, LILRB/IFITM3
 '6':'NK cell',          # NKG7, GNLY, PRF1, GZMB
 '7':'Dendritic cell',   # FCER1A, CLEC10A, HLA-DRA
 '8':'Platelet'          # PPBP, PF4, GP9
}
A.obs['cell_type']=A.obs['leiden_1.0'].map(cluster_to_type).astype(str)
assert A.obs.cell_type.notna().all()
# Put labels back in original barcode order. QC excluded 2 high-mito cells, explicitly retained as QC-failed.
raw=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/',var_names='gene_symbols')
labels=pd.Series('Low quality (QC failed)',index=raw.obs_names,dtype='object')
labels.loc[A.obs_names]=A.obs['cell_type']
out=pd.DataFrame({'barcode':raw.obs_names.astype(str),'cell_type':labels.loc[raw.obs_names].values})
out.to_csv('labels.csv',index=False)
# Save marker/annotation audit table
summary=A.obs.groupby(['leiden_1.0','cell_type'],observed=True).size().reset_index(name='n_cells')
summary.to_csv('cluster_annotation_summary.csv',index=False)
print(summary.to_string(index=False))
print('\nlabels:',out.shape,'unique barcodes:',out.barcode.nunique(),'missing:',out.isna().sum().to_dict())
print(out.cell_type.value_counts().to_string())
print('\nFile bytes',os.path.getsize('labels.csv'))