import pickle,scanpy as sc
with open('clustered.pkl','rb') as f:ad=pickle.load(f)
print(ad.shape,ad.obs.cluster.value_counts().sort_index().to_dict())
for cl in sorted(ad.obs.cluster.unique(),key=int):
 d=sc.get.rank_genes_groups_df(ad,group=cl).head(15)
 print(f"{cl}\t{sum(ad.obs.cluster==cl)}\t"+'; '.join(d.names.tolist()))