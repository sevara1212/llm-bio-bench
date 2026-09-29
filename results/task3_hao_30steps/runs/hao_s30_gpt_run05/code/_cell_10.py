import pickle,scanpy as sc,re
with open('clustered.pkl','rb') as f:ad=pickle.load(f)
for cl in [0,1,2,4,5,6,7,8,12,17,18]:
 d=sc.get.rank_genes_groups_df(ad,group=str(cl))
 g=[x for x in d.names if not re.match(r'^(RP[SL]|MT-|MALAT1|B2M|HLA-|NEAT1|TMSB)',x)][:25]
 print(cl, ', '.join(g))