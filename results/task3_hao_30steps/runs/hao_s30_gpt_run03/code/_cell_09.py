import scanpy as sc, numpy as np
ad=sc.read_h5ad('processed.h5ad'); k='leiden_0.6'
for c in ['0','4','13','18']:
 a=ad[ad.obs[k]==c]
 print('\n',c,a.n_obs)
 # raw avg (not scaled) and detection
 X=a.raw.X
 mean=np.asarray(X.mean(axis=0)).ravel(); frac=np.asarray((X>0).mean(axis=0)).ravel()
 # high expression markers excluding housekeeping
 order=np.argsort(mean)[::-1]
 out=[]
 for i in order:
  g=a.var_names[i]
  if not (g.startswith('RP') or g.startswith('MT-') or g in ['MALAT1','B2M','TMSB4X','TPT1','EEF1A1']):
   out.append(f'{g}({mean[i]:.1f},{frac[i]:.0%})')
  if len(out)==25:break
 print(', '.join(out))