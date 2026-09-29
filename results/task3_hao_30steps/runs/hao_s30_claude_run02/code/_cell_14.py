import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('adata_full_lognorm_clustered.h5ad')

df = sc.get.obs_df(adata, keys=['CD4','CD8A','CD8B','leiden'])
for cl in ['0','1','16']:
    sub = df[df.leiden==cl]
    print(cl, 'n=',len(sub))
    print(' CD4+ frac:', (sub.CD4>0.5).mean(), ' CD8A+ frac:', (sub.CD8A>0.5).mean(), ' CD8B+ frac:', (sub.CD8B>0.5).mean())