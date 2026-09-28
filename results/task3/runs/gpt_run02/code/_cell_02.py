import scipy.io as sio, pandas as pd, os
p='filtered_gene_bc_matrices/hg19'
M=sio.mmread(p+'/matrix.mtx').tocsr(); genes=pd.read_csv(p+'/genes.tsv',sep='\t',header=None); bc=pd.read_csv(p+'/barcodes.tsv',sep='\t',header=None)
print('matrix genes x cells',M.shape,'nnz',M.nnz,'total',M.sum())
print(genes.head(),genes.shape, 'unique symbols',genes[1].nunique())
print(bc.head(),bc.shape, 'barcodes unique',bc[0].nunique())
print('counts/cell', M.sum(0).min(),M.sum(0).max(),float(M.sum(0).mean()))