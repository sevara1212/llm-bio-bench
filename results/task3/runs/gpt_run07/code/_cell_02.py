import os
for root,dirs,files in os.walk('filtered_gene_bc_matrices'):
 print(root, dirs, files)