import os
for root, dirs, files in os.walk('filtered_gene_bc_matrices'):
 print(root, dirs, [(f,os.path.getsize(os.path.join(root,f))) for f in files])