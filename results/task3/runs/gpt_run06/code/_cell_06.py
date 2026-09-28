import scipy.io as sio,pandas as pd,numpy as np
p='filtered_gene_bc_matrices/hg19/'; g=np.array(pd.read_csv(p+'genes.tsv',sep='\t',header=None)[1].astype(str).tolist()); X=sio.mmread(p+'matrix.mtx').tocsr().T.tocsr()
n=np.asarray(X.sum(1)).ravel(); ng=X.getnnz(axis=1); mt=np.asarray(X[:,np.char.startswith(g,'MT-')].sum(1)).ravel()/n*100
for q in [200,300,400,500]: print(q, sum(ng>=q), 'and mt<10',sum((ng>=q)&(mt<10)))
print('cells >2500 genes',sum(ng>2500),' >8000 UMIs',sum(n>8000))