import scipy.io as sio, pandas as pd, numpy as np
from scipy import sparse
p='filtered_gene_bc_matrices/hg19/'
g=pd.read_csv(p+'genes.tsv',sep='\t',header=None)
b=pd.read_csv(p+'barcodes.tsv',sep='\t',header=None)
M=sio.mmread(p+'matrix.mtx').tocsr() # genes x cells
print('matrix',M.shape,'nnz',M.nnz,'genes',g.head().to_string(index=False),'barcodes',b.shape)
print('unique symbols',g[1].nunique(),'total',g.shape[0])
print('counts/cell',np.percentile(np.asarray(M.sum(0)).ravel(),[0,1,5,25,50,75,95,99,100]))
print('genes/cell',np.percentile(np.diff(M.tocsc().indptr),[0,1,5,25,50,75,95,99,100]))