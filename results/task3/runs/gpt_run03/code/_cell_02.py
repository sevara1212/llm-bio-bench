import os
for root,dirs,files in os.walk('filtered_gene_bc_matrices'):
 print(root, dirs, files)
 for f in files: print(' ',f,os.path.getsize(os.path.join(root,f)))