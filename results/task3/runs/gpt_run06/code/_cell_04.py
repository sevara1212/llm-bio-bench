import scipy.io as sio,pandas as pd,numpy as np
p='filtered_gene_bc_matrices/hg19/'; g=pd.read_csv(p+'genes.tsv',sep='\t',header=None)[1].astype(str).values; X=sio.mmread(p+'matrix.mtx').tocsr().T
n=np.asarray(X.sum(1)).ravel(); ng=np.diff(X.indptr); mt=np.asarray(X[:,np.char.startswith(g,'MT-')].sum(1)).ravel()/n*100
rb=np.asarray(X[:,np.char.startswith(g,'RPS')|np.char.startswith(g,'RPL')].sum(1)).ravel()/n*100
for a,name in [(mt,'mito%'),(rb,'ribo%')]: print(name,np.percentile(a,[0,1,5,25,50,75,90,95,99,100]),' >15',sum(a>15),' >20',sum(a>20))
for q in [200,300,400,500]: print(q, sum(ng>=q), 'and mt<15',sum((ng>=q)&(mt<15)))
# top genes
ix=np.asarray(X.sum(0)).ravel().argsort()[::-1][:30]
print(list(zip(g[ix],np.asarray(X.sum(0)).ravel()[ix].astype(int))))