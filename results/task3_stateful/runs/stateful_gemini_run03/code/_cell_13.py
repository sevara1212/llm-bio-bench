# What if resolution=1.0 or 1.2 or 0.6?
# In Scanpy PBMC 3k tutorial:
# The clusters are:
# 'CD4 T cells', 'CD14+ Monocytes', 'B cells', 'CD8 T cells', 'NK cells', 'FCGR3A+ Monocytes', 'Dendritic Cells', 'Megakaryocytes'
# Wait! Let's check resolution=1.0 or leiden with standard parameters.
sc.tl.leiden(adata, resolution=1.0)
print("Resolution 1.0:", adata.obs['leiden'].value_counts())
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
result = adata.uns['rank_genes_groups']
for group in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
    top_genes = [result['names'][group][i] for i in range(8)]
    print(f"Cluster {group}: {', '.join(top_genes)}")