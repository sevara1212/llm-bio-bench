import scanpy as sc
ad=sc.read_h5ad('processed_markers.h5ad'); r=ad.uns['rank_genes_groups']['names']
for g in r.dtype.names:
 vals=[]
 for x in r[g]:
  if not x.startswith(('RPS','RPL','MT-','HLA-')): vals.append(x)
 print(g, '; '.join(vals[:25]))