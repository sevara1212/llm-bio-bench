import os
for root, ds, fs in os.walk('filtered_gene_bc_matrices'):
 print(root, ds, fs)
 for f in fs: print(' ',f,os.path.getsize(os.path.join(root,f)))