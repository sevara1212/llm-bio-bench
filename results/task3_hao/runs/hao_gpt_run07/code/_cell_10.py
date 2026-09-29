import scanpy as sc, pandas as pd, numpy as np
ad=sc.read_h5ad('processed_pbmc.h5ad'); key='leiden_0.6'
labels={
'0':'Cycling T/NK cell','1':'Erythroid progenitor','2':'CD16+ monocyte','3':'CD14+ monocyte','4':'Platelet','5':'Conventional dendritic cell (cDC2)','6':'CD161+ memory T cell','7':'Erythrocyte','8':'CD4 memory T cell','9':'Memory B cell','10':'Naive B cell','11':'Naive CD4 T cell','12':'CD8 T cell','13':'Plasmacytoid dendritic cell','14':'Conventional dendritic cell (cDC1)','15':'Cytotoxic T cell','16':'Natural killer cell','17':'Gamma-delta T cell','18':'CD56bright natural killer cell','19':'Plasma cell','20':'GZMK+ memory T cell','21':'Plasmacytoid dendritic cell'}
assert set(ad.obs[key].cat.categories)==set(labels)
ad.obs['cell_type']=ad.obs[key].map(labels).astype(str)
# output all original barcodes, explicitly retaining a QC status for records outside filtered data
raw=sc.read_h5ad('raw_counts.h5ad')
out=pd.DataFrame({'barcode':raw.obs_names.astype(str)})
out['cell_type']=out.barcode.map(ad.obs['cell_type'])
# Cell types intentionally only supplied for QC-passing cells.  Flag excluded observations explicitly.
out['cell_type']=out.cell_type.fillna('QC-excluded (low complexity or putative doublet)')
out.to_csv('labels.csv',index=False)
# concise reporting table
report=(ad.obs.groupby(key,observed=True).size().rename('n_cells').to_frame())
report['cell_type']=report.index.map(labels)
report['representative_markers']=['MKI67, TYMS, TRAC, NKG7','PRSS57, GATA1, KLF1','LST1, FCGR3A, MS4A7','S100A8, S100A9, LYZ, CD14','PF4, PPBP, TUBB1','FCER1A, CD74, HLA-DRA','KLRB1, IL7R, GZMK','HBA1, HBB','IL7R, LTB, TRAC','MS4A1, CD79A, BANK1','TCL1A, IGHM, IGHD','CCR7, LTB, TCF7','CD8B, CD8A, TRAC','GZMB, GZMB, TCF4, LILRA4','CLEC9A, BATF3, CD74','CCL5, NKG7, CD3D','GNLY, PRF1, KLRD1','TRDC, GATA3, TNFRSF4','XCL1, XCL2, KLRD1','MZB1, JCHAIN, TNFRSF17','GZMK, CD27, TRAC','LILRA4, GZMB, TCF4']
report.to_csv('cluster_annotation_report.csv',index_label='cluster')
print(report.to_string())
print('\nlabels shape:',out.shape,'unique barcodes:',out.barcode.nunique(),'missing:',out.cell_type.isna().sum())
print(out.cell_type.value_counts().to_string())
print('\nQC-excluded:',(out.cell_type.str.startswith('QC-excluded')).sum())