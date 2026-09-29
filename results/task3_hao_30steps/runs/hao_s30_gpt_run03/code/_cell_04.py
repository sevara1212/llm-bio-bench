import pandas as pd
x=pd.read_csv('qc_metrics.csv',index_col=0)
for mt in [10,15,20,25,30,40]:
 m=(x.n_genes_by_counts>=500)&(x.pct_counts_mt<mt)&(x.pct_counts_hb<20)&(x.n_genes_by_counts<5500)
 print(mt,m.sum(), 'mtfail',(x.pct_counts_mt>=mt).sum(),'hbfail',(x.pct_counts_hb>=20).sum(),'highgene',(x.n_genes_by_counts>=5500).sum())