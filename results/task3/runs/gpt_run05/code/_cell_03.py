import scipy.io as sio, pandas as pd, numpy as np
from scipy import sparse
p='filtered_gene_bc_matrices/hg19/'
genes=pd.read_csv(p+'genes.tsv',sep='\t',header=None)
bc=pd.read_csv(p+'barcodes.tsv',sep='t',header=None)[0].astype(str).values
X=sio.mmread(p+'matrix.mtx').tocsr() # genes x cells
print('shape genes x cells',X.shape,'nnz',X.nnz)
print(genes.head(), genes.tail(), sep='\n')
print('barcodes',bc[:3],bc[-3:], 'unique',len(set(bc)))
cs=np.asarray(X.sum(0)).ravel(); gs=np.asarray((X>0).sum(0)).ravel()
for a in [cs,gs]: print(np.percentile(a,[0,1,5,25,50,75,95,99,100]), 'mean',a.mean())