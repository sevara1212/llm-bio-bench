import pickle, pandas as pd, numpy as np
with open('clustered.pkl','rb') as f: ad=pickle.load(f)
genes=['CD3D','CD3E','TRAC','IL7R','CCR7','LTB','MALAT1','NKG7','GNLY','MS4A1','CD79A','CD74','LYZ','S100A8','FCGR3A','FCER1A','CLEC9A','GZMB','TRDC','MKI67','HBB','PRSS57','GATA2','KIT','CD34','JCHAIN','MZB1','LILRA4']
raw=ad.raw.to_adata();x=raw[:,genes].X
if hasattr(x,'toarray'):x=x.toarray()
df=pd.DataFrame(x,index=ad.obs.cluster,columns=genes).groupby(level=0,observed=False).mean()
print(df.round(2).to_string())